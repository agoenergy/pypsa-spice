// SPDX-FileCopyrightText: PyPSA-SPICE Developers
// SPDX-License-Identifier: GPL-2.0-or-later

import { useEffect, useState } from "react";
import { createModelRunMonitor } from "./modelRunMonitor";
import type { ModelRunMonitor, ModelRunMonitorState } from "./types";

export default function useModelRun(): ModelRunMonitor {
  const [state, setState] = useState<ModelRunMonitorState>({ run: null, error: "", loading: true });
  const [monitor] = useState(() => createModelRunMonitor(setState));

  useEffect(() => {
    monitor.start();
    return () => monitor.stop();
  }, [monitor]);

  return { ...state, updateRun: monitor.updateRun };
}
