# CircosContacts 1.0 Ticket Set

## CC-101: State Model v2 + Editable Plot Title

### Problem
The plot lacks an editable title that persists across exports and sessions.

### Scope
- Add `plotTitle` to runtime state.
- Add title text input in left panel.
- Render title centered above plot.
- Include title in SVG export.
- Include title in session save/load.
- Introduce `session_schema: v2` while preserving load support for v1 session files.

### Acceptance Criteria
- Editing title updates plot immediately.
- `Download SVG` contains same title text.
- New session JSON includes title and reloads it correctly.
- Old v1 session JSON still loads without errors.

### Non-goals
- Advanced typography/styling controls for title.

---

## CC-102: Auto-fit Layout (Remove Manual Zoom)

### Problem
Manual zoom introduces clipping and extra control complexity.

### Scope
- Remove zoom controls.
- Auto-fit plot to viewport/canvas area.
- Keep plot centered by default.
- Ensure tooltip and menus remain visible within viewport.

### Acceptance Criteria
- No zoom controls shown.
- Plot is fully visible on initial load at common desktop sizes.
- Tooltip/menu never render outside visible panel bounds.

### Non-goals
- User-controlled zoom restoration.

---

## CC-103: Reproducible Save HTML

### Problem
`Save HTML` currently does not guarantee full-state reproducibility without extra steps.

### Scope
- Save current runtime state payload directly in HTML snapshot.
- Reopening saved HTML restores exact state.
- Keep `Save Session` behavior for transferring settings to other datasets.

### Acceptance Criteria
- Saved HTML reopens with matching controls, chain settings, title, selections, callouts, and colors.
- `Save Session` remains available and unchanged in purpose.

### Non-goals
- Multi-file export bundle format.

---

## CC-104: Arc-following Labels + Selection-driven Link Coloring

### Problem
Chain labels are not path-aligned; links cannot be recolored by selected regions.

### Scope
- Render chain labels on arc-following text paths, offset from arcs.
- Orient labels for readability based on hemisphere.
- Add selection option to recolor contact links for links intersecting selection residues.
- Persist link color overrides in session and SVG export.

### Acceptance Criteria
- Labels follow arc contour and remain readable top/bottom.
- Link recolor overrides apply only to intended links.
- SVG export matches on-screen label orientation and link colors.

### Non-goals
- Per-link manual color editing.

---

## CC-105: Stability/UX Hardening

### Problem
Editing workflows need stronger recoverability and provenance.

### Scope
- Add undo/redo for selection and annotation actions.
- Add compact legend/state stamp (title + key settings) for export provenance.
- Improve render performance for dense datasets (debounce/cache/hit-test optimizations).
- Update docs for all 1.0 behavior changes.

### Acceptance Criteria
- Undo/redo restores previous states for key edit actions.
- Exported SVG can include legend/state stamp.
- Dense plots remain responsive compared with prior behavior.
- Documentation updated and internally consistent.

### Non-goals
- New analysis algorithms.
