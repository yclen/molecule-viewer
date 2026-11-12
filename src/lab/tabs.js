// The three app tabs: the molecule viewer, the live pulse experiment and the intensity sweep. The address's #hash
// picks one, so a link can open straight onto a tab. The two lab pages are built the first time they're shown.

import './lab.css';

const PAGES = { viewer: 'Molecule Viewer', pulse: 'Pulse experiment', sweep: 'Intensity sweep' };
const mounted = {};

async function show(page, onChange) {
  if (!(page in PAGES)) page = 'viewer';
  document.body.dataset.page = page;
  document.title = page === 'viewer' ? 'Molecule Viewer' : `${PAGES[page]} · Molecule Viewer`;
  for (const a of document.querySelectorAll('#apptabs a[data-page]')) {
    if (a.dataset.page === page) a.setAttribute('aria-current', 'page');
    else a.removeAttribute('aria-current');
  }
  for (const name of ['pulse', 'sweep']) {
    const root = document.getElementById(`page-${name}`);
    root.hidden = name !== page;
    if (mounted[name]) (await mounted[name]).shown?.(name === page);
  }
  if (page !== 'viewer' && !mounted[page]) {
    const root = document.getElementById(`page-${page}`);
    mounted[page] = (page === 'pulse' ? import('./pulse.js') : import('./sweep.js')).then((m) => m.mount(root));
    (await mounted[page]).shown?.(true);
  }
  onChange(page);
}

export function initTabs(onChange) {
  const current = () => location.hash.slice(1);
  window.addEventListener('hashchange', () => show(current(), onChange));
  show(current(), onChange);
}
