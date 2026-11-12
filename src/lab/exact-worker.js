// The exact solves, off the page's own thread so it stays responsive. A solve can't be interrupted, so to drop one the
// page ends this worker and starts a fresh one.
//
// { type: 'pulseExact', key, p, R }       -> { type: 'progress', key, share, doing } ..., then { type: 'done', key, result }
// { type: 'totals', key, p, R, rules }    -> { type: 'totals', key, counts }
// Either can answer { type: 'failed', key, error } instead.

import { pulseExact, pulseTotals } from '../sim/exact.js';

onmessage = ({ data: m }) => {
  try {
    if (m.type === 'pulseExact') {
      let last = 0;
      const progress = (share, doing) => {
        if (performance.now() - last > 100) {
          last = performance.now();
          postMessage({ type: 'progress', key: m.key, share, doing });
        }
      };
      postMessage({ type: 'done', key: m.key, result: pulseExact(m.R, m.p, m.p.window, progress) });
    } else if (m.type === 'totals') {
      const counts = pulseTotals(m.R, m.p, m.rules);
      postMessage({ type: 'totals', key: m.key, counts }, [counts.buffer]);
    }
  } catch (error) {
    postMessage({ type: 'failed', key: m.key, error: String(error.message ?? error) });
  }
};
