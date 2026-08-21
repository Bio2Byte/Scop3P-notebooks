"""Waiting-time handling shared by every protocol.

These protocols do slow things: a Bio2Byte prediction takes ~16s, a RIN build runs a
KD-tree over a few hundred residues, TM-align shells out, and every annotation source is a
network round trip. Without explicit handling the user gets a page that looks broken --
nothing moves, so they click again, and the second click queues another run of the same
work.

The mechanics here are not obvious, and were established by measuring a running app rather
than read from documentation, so they are written down.

**A synchronous handler cannot show a busy state at all.** Everything an effect writes --
a status message, a disabled button -- is queued and flushed only after the effect returns.
Measured: a 3-second handler that sets "Working..." and disables its button shows *neither*
while it runs; the page sits unchanged and then jumps straight to the finished state. So
the work has to leave the event loop.

**An async function is not enough either.** ``reactive.extended_task`` runs its coroutine
on the event loop, so a blocking call inside ``async def`` still pins it. The blocking work
must go to a thread with ``asyncio.to_thread``. That is also what stops one slow upstream
freezing *every* connected session -- a real fault seen here when the EBI Proteins API
accepted a connection and then never answered.

**The button's busy state needs driving.** ``input_task_button`` disables itself client-side
the instant it is clicked, which is what actually prevents a double click. With
``auto_reset=True`` it re-enables on the next flush -- measured at 177ms into a 3-second
task, so the button is live again while the work is still running. These buttons therefore
use ``auto_reset=False`` and are reset explicitly when the task ends.

Timeline of the resulting behaviour, measured:

===========  ==========================================================
~2ms         button disabled (client-side, so a second click is refused)
~250ms       label becomes "Working...", status card shows what is running
end of work  status shows the outcome, button re-enabled -- on success *and* on error
===========  ==========================================================
"""

from __future__ import annotations

import asyncio
from typing import Any, Callable

from shiny import reactive, ui

#: Label a task button wears while its work runs.
BUSY_LABEL = "Working…"

# Both button helpers forward **kwargs to ui.input_task_button. That is load-bearing:
# several buttons are rendered inside a @render.ui output and pass `disabled=` to reflect
# whether their prerequisite has been met -- Run TM-align until structures are loaded,
# Compare / Align until both networks exist. Swallowing that keyword raised a TypeError
# *at render time*, which replaced the whole output with an error message and made the
# button disappear. Construction-time smoke tests never render those outputs, so nothing
# caught it; test_busy_state now checks the signature against every call site's keywords.


def busy_indicators():
    """Page-level spinners for outputs that are recalculating.

    Include once per app, inside the page. Unlike a status message this *does* work with
    synchronous handlers, because Shiny marks dependent outputs busy in the browser
    without waiting for the server.
    """
    return ui.busy_indicators.use(spinners=True, pulse=True)


def task_button(
    input_id: str,
    label: str,
    *,
    class_: str = "",
    label_busy: str = BUSY_LABEL,
    **kwargs: Any,
):
    """A button that guards itself while an ordinary synchronous handler runs.

    Measured against a 3-second blocking handler: disabled 3ms after the click (so a
    second click 500ms later is refused), spinner and busy label showing from ~760ms, and
    restored when the handler returns. Three of the four things wanted, from a one-line
    change, with no restructuring.

    What it cannot do is show a *live* status message: that is queued until the handler
    returns, because the handler holds the event loop. Where the wait is long enough that
    the user needs to be told what is happening, use ``background_task_button`` and move
    the work off the loop.

    ``auto_reset`` is left at its default, so the button restores itself on the flush that
    follows the handler. No ``finish_task`` call is needed or wanted.
    """
    return ui.input_task_button(
        input_id,
        label,
        label_busy=label_busy,
        class_=class_,
        **kwargs,
    )


def background_task_button(
    input_id: str,
    label: str,
    *,
    class_: str = "",
    label_busy: str = BUSY_LABEL,
    **kwargs: Any,
):
    """A button for work that runs off the event loop, held busy for its whole duration.

    ``auto_reset=False`` is the point: with the default the button re-enables on the next
    flush -- measured at 177ms into a 3-second task -- so it would be live again while the
    work was still running.

    The contract that comes with that: **every exit path must call ``finish_task``**,
    including early returns like "set an accession first". A missed path leaves the button
    dead for the rest of the session with nothing on screen to explain it. A test walks the
    source to check each of these buttons has a matching ``finish_task``.
    """
    return ui.input_task_button(
        input_id,
        label,
        label_busy=label_busy,
        auto_reset=False,
        class_=class_,
        **kwargs,
    )


#: Class that makes a control visibly inert. Defined here so every protocol gates the
#: same way and the CSS can live in one stylesheet.
INERT_CLASS = "scop3p-inert"

INERT_CSS = """
.scop3p-inert {
  opacity: 0.45;
  pointer-events: none;
  user-select: none;
}
"""


def gate(control, *, ready: bool, hint: str = ""):
    """Show a control as unavailable until its prerequisite is met.

    A plain action button can be gated with ``disabled=True``. A **task button cannot**:
    ``bslib-task-button`` owns the ``disabled`` property as part of its ready/busy state
    machine and drops the attribute when it initialises, so the server's flag is silently
    discarded and the button renders fully clickable. Verified in the browser -- the
    rendered markup carries no ``disabled`` attribute at all.

    So the control is wrapped instead: ``pointer-events: none`` refuses the click and the
    reduced opacity shows why. An optional hint says what has to happen first, because a
    greyed control with no explanation is its own kind of dead end.
    """
    if ready:
        return control
    children = [control]
    if hint:
        children.append(ui.p(hint, class_="scop3p-note scop3p-inert-hint"))
    return ui.div(*children, class_=INERT_CLASS)


def background(function: Callable[..., Any]):
    """Wrap a blocking callable as an extended task that runs off the event loop.

    Call inside ``server()`` -- an extended task belongs to one session. The wrapped
    callable must not touch reactive values or ``input``: it runs in a worker thread,
    outside any reactive context. Pass what it needs as plain arguments and apply the
    result in an effect that reads ``.result()``.
    """

    @reactive.extended_task
    async def _task(*args: Any, **kwargs: Any) -> Any:
        return await asyncio.to_thread(function, *args, **kwargs)

    return _task


def finish_task(input_id: str) -> None:
    """Re-enable a task button.

    Must be called on *every* exit path. The click already disabled the button in the
    browser, so a handler that returns early -- "set an accession first" -- without this
    leaves the control dead for the rest of the session with nothing to explain why.
    """
    ui.update_task_button(input_id, state="ready")


def task_outcome(task, *, on_success, on_error, on_finished=None) -> bool:
    """Apply a finished task's result, once, in a collector effect.

    Returns True when the task had ended (either way), which is when the button should be
    released. A task that is still running, or has never been started, does nothing --
    reading ``.result()`` in that state raises, and an effect that raises would abort the
    flush that is carrying the "Working..." message to the browser.
    """
    state = task.status()
    if state == "success":
        on_success(task.result())
    elif state == "error":
        try:
            task.result()
        except Exception as error:  # noqa: BLE001 - reported, not swallowed
            on_error(error)
    else:
        return False
    if on_finished is not None:
        on_finished()
    return True
