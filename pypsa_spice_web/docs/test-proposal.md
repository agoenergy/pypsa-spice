# Frontend test proposal

Written 2026-09-12, for discussion. Nothing here is implemented yet.

This picks up [Jia General 5](ui-review-jia.md#general), which settled on tests grouped
by behaviour, alongside the code that owns that behaviour, plus dedicated component
tests where a component has meaningful interactions. The reply there also noted that
broader behaviour coverage remains open. This is a concrete first slice of it.

## What prompted it

The simplify pass on 2026-09-12 removed the responsive stylesheets and the mobile
navigation they drove, untangled the `ChartSurface` override tail, and pulled four
shared modules out of duplicated code:

| Module | Replaces |
| --- | --- |
| `src/useDismiss.ts` | The Escape handler written out five times |
| `src/useNearViewport.ts` | New. Gates chart requests behind an `IntersectionObserver` |
| `src/components/ResultsToc.tsx` | Three copies of the same jump-to popover |
| `src/components/ChartFeedback.tsx` | Loading, error, and heading blocks duplicated across both chart cards |

Three of them own behaviour that fails without throwing, which is the case for
dedicated tests. None of them has a test today.

## Why nothing covers them now

The suite is eleven `.ts` files and 44 tests, all pure logic. There is no jsdom, no
`@testing-library/react`, and no `environment` setting in `vite.config.ts`. React hooks
are not reachable from that setup, so the gap is structural rather than an oversight.

## Proposed tests

Four files, roughly 18 tests, each beside the code it covers.

| File | Covers | Tests |
| --- | --- | --- |
| `src/useNearViewport.test.ts` | Lazy chart loading gate | 5 |
| `src/useDismiss.test.ts` | Escape and outside-click dismissal | 5 |
| `src/components/ChartFeedback.test.ts` | `useDelayedFlag` spinner timing | 4 |
| `src/components/ResultsToc.test.tsx` | Popover open, close, and empty guard | 4 |

### `useNearViewport`

The one I care most about. Starts false, latches true on first intersection, stays true
after the element scrolls away, disconnects the observer on unmount, and falls back to
true when `IntersectionObserver` is missing.

Both failure modes are silent. Stuck false means charts never load. Stuck true means the
Power page quietly returns to firing all 16 requests on mount, which is what this module
exists to prevent. Neither throws, and neither shows up in a diff review.

### `useDismiss`

Four call sites depend on it, and one gates on state: `NewScenarioDialog` ignores Escape
while a submit is in flight. Escape fires while active, does nothing while inactive, the
listener detaches on unmount, an outside pointerdown dismisses, an inside one does not.

`DataDialog` keeps its own handler, because Escape there closes open table menus first
and only then the dialog. That stays out of the shared hook and out of these tests.

### `useDelayedFlag`

The guarantee is that a spinner does not flash on fast responses. That is pure timing:
false before the delay, true after it, back to false when loading clears, and no late
flip if loading clears before the timer fires. Fake timers, so no real waiting.

### `ResultsToc`

The debatable one, and worth arguing about. It sits close to the visual wrapper category
that Jia said is not worth a test file.

The case for covering it: one component now serves Results, Scenario differences, and
Configure, so a regression hits three pages at once. It also carries real logic. It
returns null when `entries` is empty, `aria-expanded` tracks panel state, and selecting
an entry closes the panel. Tests would cover those four things and nothing else.

Happy to drop this file if Jia reads it as a wrapper.

## What I would not test

`ChartLoadingState`, `ChartErrorState`, `ChartRefreshOverlay`, and `ChartContextHeading`
in `ChartFeedback.tsx`. They are markup with no branching. Testing them would assert that
JSX renders JSX, which is the maintenance-without-value case from General 5.

## Cost

Two development dependencies, `jsdom` and `@testing-library/react` v16 or later, which is
the first version supporting React 19. The project is on React 19.2.7.

The existing suite does not have to move. Vitest 4.1.10 supports a per-file
`// @vitest-environment jsdom` docblock, verified in the installed package, so the eleven
current files keep running in the node environment at their current speed and only the new
files pay for a DOM.

## The decision underneath

These four files matter less than the precedent they set. Approving them means the project
accepts a DOM test environment and component-level tests as a normal pattern. That opens
coverage for behaviour that is untestable today and has actually broken before: dirty-state
save and discard, run polling, the deep-link race that drops `run` and `compare`, and
shared axis ranges in comparison mode.

That is the question worth the discussion. The four files are a small first instance of it.

## Recommendation

Take it, `ResultsToc` included. The hooks encode invariants that are invisible in review
and silent when they break, and the marginal cost after the first jsdom file is close to
zero.

If the answer is no DOM at all, then the honest outcome is that these four modules stay
untested. Reshaping them into pure functions to fit the current setup would produce tests
that pass while the lifecycle behaviour, which is where the risk actually sits, goes
unchecked. I would rather have no test than that one.

## Open questions

1. Is `ResultsToc` worth a test file, or does it fall under the visual wrapper exclusion?
2. Per-file `@vitest-environment` docblocks, or one global jsdom environment for the whole
   suite? Per-file keeps the current tests fast; global is less to remember.
3. Does `@testing-library/react` fit, or is there a lighter approach you would prefer for
   hook lifecycle tests?
4. If this pattern is accepted, which of the untested behaviours above comes next?
