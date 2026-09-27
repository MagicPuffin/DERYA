# Task Delegation Workflow

Read this before starting any task from `backlog.md`.

## What makes a task safe to run in one sitting
1. **One function, one file, or one well-bounded edit.**
2. **Interface already fixed** in the task description — not invented mid-task.
3. **A failing test already exists, or is itself the task.** Success =
   "this specific test goes green," never "make it work."
4. **Scope is explicit** — the task states which files may change.
   `git diff --stat` after is the check; anything outside the list is a
   reject-and-redo, not a judgment call. An exception the agent flags and
   justifies (e.g. a shared test-CMakeLists edit needed to link a new module)
   is fine — an exception it doesn't flag is not.

## Stop condition (put this in every task prompt)
> If the test isn't green after ~2-3 genuine attempts, stop and report what
> was tried and where it diverges. Do not keep rewriting.

This applies symmetrically: if the *test's expected value* looks wrong
against the spec/equation rather than the implementation being wrong, stop
and report that too — do not "fix" a test to match the code just to get to
green. (This happened for real on Task 1 — the human-supplied expected value
was arithmetically wrong, not the code. Correctly flagged instead of silently
"fixed.")

## Session sizing
One task per session, review the diff, then decide whether to queue the
next. Batch only once a module has a track record of clean single-task
sessions, and only independent tasks that don't touch shared files.

## Review checklist (run every time)
- [ ] The specified test, and only its target, went red → green
- [ ] `git diff --stat` matches "files allowed to change" (flagged exceptions
      reviewed on their own merits)
- [ ] Full test suite still passes
- [ ] No stray TODOs/stub shortcuts left in the diff
- [ ] `backlog.md` checkboxes updated for the task
- [ ] `README.md` updated if the task changed what the code can do (a new
      capability under "What it does today") or passed a milestone gate
      (the "Status" section). Detailed progress stays in `backlog.md` and
      `roadmap.md`; the README only summarizes and links to them
- [ ] One commit, message references the task
