// Saving runs. A run is a set of files (results.csv, params.json, ...), the same ones sim/ writes into a run folder.
// It's saved one of two ways:
//  - downloaded as one .zip (any browser)
//  - written into a folder the user picked once (Chrome and Edge: the File System Access API). Picking a folder that
//    Google Drive, Dropbox or OneDrive syncs puts the runs there.
// Sweep runs are also kept in this browser (IndexedDB), so the Plot list survives a reload.

const DB_NAME = 'molecule-viewer-lab';
let dbPromise = null;

function db() {
  dbPromise ??= new Promise((resolve, reject) => {
    const open = indexedDB.open(DB_NAME, 1);
    open.onupgradeneeded = () => {
      open.result.createObjectStore('runs', { keyPath: 'id' });
      open.result.createObjectStore('settings');
    };
    open.onsuccess = () => resolve(open.result);
    open.onerror = () => reject(open.error);
  });
  return dbPromise;
}

async function request(store, mode, act) {
  const tx = (await db()).transaction(store, mode);
  const req = act(tx.objectStore(store));
  return new Promise((resolve, reject) => {
    tx.oncomplete = () => resolve(req?.result);
    tx.onerror = () => reject(tx.error);
  });
}

export const listRuns = () => request('runs', 'readonly', (s) => s.getAll()).catch(() => []);
export const putRun = (run) => request('runs', 'readwrite', (s) => s.put(run)).catch(() => {});
export const deleteRun = (id) => request('runs', 'readwrite', (s) => s.delete(id)).catch(() => {});
const getSetting = (key) => request('settings', 'readonly', (s) => s.get(key)).catch(() => undefined);
const setSetting = (key, value) => request('settings', 'readwrite', (s) => s.put(value, key)).catch(() => {});

// ---- a .zip with no compression ----

const CRC = new Uint32Array(256).map((_, n) => {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  return c;
});
const crc32 = (bytes) => {
  let c = 0xffffffff;
  for (const b of bytes) c = CRC[(c ^ b) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
};

/** files: [{ name, data: string | Uint8Array }] -> a Blob holding a .zip of them. */
export function zip(files) {
  const enc = new TextEncoder();
  const now = new Date();
  const time = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  const date = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate();
  const parts = [];
  const central = [];
  let offset = 0;
  for (const file of files) {
    const name = enc.encode(file.name);
    const data = typeof file.data === 'string' ? enc.encode(file.data) : file.data;
    const crc = crc32(data);
    const local = new DataView(new ArrayBuffer(30));
    [[0, 0x04034b50, 4], [4, 20, 2], [6, 0x0800, 2], [8, 0, 2], [10, time, 2], [12, date, 2], [14, crc, 4],
      [18, data.length, 4], [22, data.length, 4], [26, name.length, 2], [28, 0, 2]]
      .forEach(([at, v, size]) => (size === 4 ? local.setUint32(at, v, true) : local.setUint16(at, v, true)));
    const entry = new DataView(new ArrayBuffer(46));
    [[0, 0x02014b50, 4], [4, 20, 2], [6, 20, 2], [8, 0x0800, 2], [10, 0, 2], [12, time, 2], [14, date, 2], [16, crc, 4],
      [20, data.length, 4], [24, data.length, 4], [28, name.length, 2], [30, 0, 2], [32, 0, 2], [34, 0, 2], [36, 0, 2],
      [38, 0, 4], [42, offset, 4]]
      .forEach(([at, v, size]) => (size === 4 ? entry.setUint32(at, v, true) : entry.setUint16(at, v, true)));
    parts.push(local, name, data);
    central.push(entry, name);
    offset += 30 + name.length + data.length;
  }
  const size = central.reduce((s, p) => s + p.byteLength, 0);
  const end = new DataView(new ArrayBuffer(22));
  [[0, 0x06054b50, 4], [8, files.length, 2], [10, files.length, 2], [12, size, 4], [16, offset, 4]]
    .forEach(([at, v, s]) => (s === 4 ? end.setUint32(at, v, true) : end.setUint16(at, v, true)));
  return new Blob([...parts, ...central, end], { type: 'application/zip' });
}

export function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 10_000);
}

// ---- where saved runs go ----

export const canPickFolder = typeof window !== 'undefined' && 'showDirectoryPicker' in window;
let folder = null; // a FileSystemDirectoryHandle, once picked
const listeners = new Set();
const announce = () => listeners.forEach((f) => f());
getSetting('folder').then((handle) => {
  if (handle) {
    folder = handle;
    announce();
  }
});

export async function pickFolder() {
  const handle = await window.showDirectoryPicker({ id: 'tta-runs', mode: 'readwrite' });
  folder = handle;
  await setSetting('folder', handle);
  announce();
}

export async function useDownloads() {
  folder = null;
  await setSetting('folder', null);
  announce();
}

export const destination = () => (folder ? folder.name : null);
export const onFolderChange = (f) => listeners.add(f);

/**
 * Whether the page may read and write the picked folder. The browser remembers the folder between visits but asks
 * again for permission; `ask` shows that prompt, which only works straight after a click.
 */
export async function folderAccess(ask = false) {
  if (!folder) return false;
  try {
    if ((await folder.queryPermission({ mode: 'readwrite' })) === 'granted') return true;
    return ask && (await folder.requestPermission({ mode: 'readwrite' })) === 'granted';
  } catch {
    return false;
  }
}

/**
 * A new run folder `name` inside the picked folder (`name_2`, `name_3`... if taken, as sim/ does). Returns
 * { name, where, write(fileName, data) }; write replaces the file each time, so a run can be saved as it goes.
 */
export async function makeRunFolder(name) {
  let dir = name;
  for (let k = 2; ; k++) {
    try {
      await folder.getDirectoryHandle(dir);
      dir = `${name}_${k}`;
    } catch {
      break;
    }
  }
  const sub = await folder.getDirectoryHandle(dir, { create: true });
  return {
    name: dir,
    where: `${folder.name}/${dir}`,
    async write(fileName, data) {
      const writable = await (await sub.getFileHandle(fileName, { create: true })).createWritable();
      await writable.write(data);
      await writable.close();
    },
    /** Remove this run folder again: only for one this page just made and never filled. */
    remove: () => folder.removeEntry(dir, { recursive: true }),
  };
}

const readText = async (dir, name) => (await (await dir.getFileHandle(name)).getFile()).text();

/**
 * Every intensity sweep in the picked folder: each subfolder whose params.json says mode "intensity" and that has a
 * results.csv. The website's sweeps and the Python app's (sim/data) have the same layout, so both show.
 */
export async function listFolderRuns() {
  if (!(await folderAccess())) return [];
  const runs = [];
  for await (const [name, handle] of folder.entries()) {
    if (handle.kind !== 'directory') continue;
    try {
      const params = JSON.parse(await readText(handle, 'params.json'));
      if (params.mode !== 'intensity') continue;
      const rows = fromCsv(await readText(handle, 'results.csv')).map((r) => ({ ...r, method: String(r.method ?? 'exact') }));
      runs.push({ id: name, date: params.date ?? '', params, rows, where: 'folder' });
    } catch {
      // not a run folder, or one still being written
    }
  }
  return runs;
}

/** Save a run's files: into the picked folder (as a subfolder `name`), or as `name`.zip. Returns where it went. */
export async function saveRun(name, files) {
  if (folder) {
    try {
      if (!(await folderAccess(true))) throw new Error('no permission');
      const run = await makeRunFolder(name);
      for (const file of files) await run.write(file.name, file.data);
      return run.where;
    } catch (error) {
      if (error.name === 'AbortError') throw error;
      // the folder went away or permission was refused: fall back to a download
    }
  }
  download(`${name}.zip`, zip(files));
  return `${name}.zip in your downloads`;
}

/** A small "Saves go to …" control, kept in step on every page. */
export function saveLocation() {
  const box = document.createElement('div');
  box.className = 'save-location';
  const render = () => {
    const where = destination();
    box.innerHTML = `<span>Saves go to</span> <b>${where ? `📁 ${escapeHtml(where)}` : 'Downloads (.zip)'}</b>`;
    if (canPickFolder) {
      const link = (text, title, act) => {
        const b = document.createElement('button');
        b.className = 'link';
        b.textContent = text;
        b.title = title;
        b.onclick = act;
        box.append(b);
      };
      const choose = () => pickFolder().catch(() => {});
      if (where) {
        link('Change…', 'Pick another folder', choose);
        link('Use downloads', 'Download each run as a .zip instead', useDownloads);
      } else {
        link('Choose a folder…', 'Save runs straight into a folder on this computer, and open the sweeps already in it. '
          + 'A folder Google Drive or Dropbox syncs works too (Chrome and Edge)', choose);
      }
    }
  };
  listeners.add(render);
  render();
  return box;
}

export const escapeHtml = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);

/** Rows of objects -> CSV text, with the columns of the first row. */
export function toCsv(rows) {
  if (!rows.length) return '';
  const keys = Object.keys(rows[0]);
  const cell = (v) => {
    if (v == null || (typeof v === 'number' && Number.isNaN(v))) return 'nan';
    const s = typeof v === 'number' ? String(v) : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return `${keys.join(',')}\n${rows.map((r) => keys.map((k) => cell(r[k])).join(',')).join('\n')}\n`;
}

/** CSV text -> rows of objects, numbers as numbers. */
export function fromCsv(text) {
  const lines = text.replace(/\r/g, '').split('\n').filter(Boolean);
  const split = (line) => {
    const out = [];
    let cur = '';
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (quoted) {
        if (c === '"' && line[i + 1] === '"') cur += line[i++];
        else if (c === '"') quoted = false;
        else cur += c;
      } else if (c === '"') quoted = true;
      else if (c === ',') {
        out.push(cur);
        cur = '';
      } else cur += c;
    }
    out.push(cur);
    return out;
  };
  const keys = split(lines[0]);
  return lines.slice(1).map((line) => {
    const cells = split(line);
    return Object.fromEntries(keys.map((k, i) => {
      const v = cells[i] ?? '';
      const num = v === '' || v === 'nan' ? NaN : Number(v);
      return [k, Number.isNaN(num) && v !== '' && v !== 'nan' ? v : num];
    }));
  });
}

export const stamp = (d = new Date()) => {
  const pad = (n) => String(n).padStart(2, '0');
  return { day: `${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())}`,
    time: `${pad(d.getHours())}:${pad(d.getMinutes())}:${pad(d.getSeconds())}` };
};
