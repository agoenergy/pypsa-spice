const dirtyEditors = new Set<string>();

export function setEditorDirty(id: string, dirty: boolean) {
  if (dirty) dirtyEditors.add(id);
  else dirtyEditors.delete(id);
}

export function confirmDiscardChanges() {
  return dirtyEditors.size === 0 || window.confirm("Discard unsaved changes before changing workspace context?");
}
