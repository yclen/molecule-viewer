import * as THREE from 'three';
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { getElement } from './elements.js';
import { DEPTH_LAYER, GlowLayer } from './glow.js';

/** Display styles. Radii are in Å. */
export const STYLES = {
  ballstick: { label: 'Ball & stick', atomRadius: (el) => el.vdw * 0.25, bondRadius: 0.1 },
  sticks: { label: 'Sticks', atomRadius: () => 0.15, bondRadius: 0.15 },
  spacefill: { label: 'Space-filling', atomRadius: (el) => el.vdw, bondRadius: 0 },
  wireframe: { label: 'Wireframe', atomRadius: () => 0.07, bondRadius: 0, lines: true },
};

const Y_AXIS = new THREE.Vector3(0, 1, 0);

export class MoleculeViewer {
  constructor(container) {
    this.container = container;
    this.style = 'ballstick';
    this.showHydrogens = true;
    this.atoms = [];
    this.bonds = new Uint32Array(0);
    this.visible = []; // instance index -> atom index
    this.radii = new Float32Array(0);
    this.center = new THREE.Vector3();
    this.meshes = [];
    this.atomMesh = null;
    this.hovered = -1;
    this.onHover = null; // (atomIndex, atom, clientX, clientY) => void

    // --- renderer & scene -------------------------------------------------
    const renderer = (this.renderer = new THREE.WebGLRenderer({
      antialias: true,
      alpha: true,
      preserveDrawingBuffer: true, // lets us grab screenshots later
    }));
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);

    this.scene = new THREE.Scene();
    const pmrem = new THREE.PMREMGenerator(renderer);
    this.scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    this.scene.environmentIntensity = 0.55;
    pmrem.dispose();

    this.camera = new THREE.PerspectiveCamera(35, 1, 0.1, 1000);
    // Key light rides along with the camera so the molecule is always lit from the viewer's upper-left.
    const key = new THREE.DirectionalLight(0xffffff, 2.0);
    key.position.set(-1.2, 1.5, 2.5);
    key.target.position.set(0, 0, -1);
    this.camera.add(key, key.target);
    this.scene.add(this.camera);
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.25));

    // Trackball = free rotation in every direction (no "up" lock), like Avogadro.
    const controls = (this.controls = new TrackballControls(this.camera, renderer.domElement));
    controls.rotateSpeed = 3.5;
    controls.zoomSpeed = 1.4;
    controls.panSpeed = 0.6;
    controls.dynamicDampingFactor = 0.15;

    // pivot = principal-axis rotation, group = translation to put the centroid at the origin
    this.pivot = new THREE.Group();
    this.group = new THREE.Group();
    this.pivot.add(this.group);
    this.scene.add(this.pivot);

    // --- shared geometry & materials -------------------------------------
    this.sphereGeos = {};
    this.cylGeo = new THREE.CylinderGeometry(1, 1, 1, 20, 1, true);
    this.atomMat = new THREE.MeshStandardMaterial({ roughness: 0.32, metalness: 0.0 });
    this.bondMat = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.0 });
    this.lineMat = new LineMaterial({ linewidth: 2.5, vertexColors: true, worldUnits: false });

    this.halo = new THREE.Mesh(
      this._sphereGeo(0),
      new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.28, depthWrite: false }),
    );
    this.halo.visible = false;
    this.halo.renderOrder = 1;
    this.group.add(this.halo);

    // excited-state glow (drawn in its own passes, see glow.js)
    this.glow = new GlowLayer(renderer);
    this.group.add(this.glow.mesh);
    this.clock = new THREE.Clock();

    // --- hover picking ----------------------------------------------------
    this.raycaster = new THREE.Raycaster();
    this.pointer = { x: 0, y: 0, inside: false, dirty: false, down: false };
    const canvas = renderer.domElement;
    canvas.addEventListener('pointermove', (e) => {
      Object.assign(this.pointer, { x: e.clientX, y: e.clientY, inside: true, dirty: true });
    });
    canvas.addEventListener('pointerleave', () => {
      this.pointer.inside = false;
      this._setHover(-1);
    });
    canvas.addEventListener('pointerdown', (e) => {
      this.pointer.down = true;
      this.pointer.downAt = e.button === 0 ? { x: e.clientX, y: e.clientY } : null;
      this._setHover(-1);
    });
    window.addEventListener('pointerup', (e) => {
      this.pointer.down = false;
      this.pointer.dirty = true;
      // A click (not a rotate-drag) on the canvas picks an atom for measuring.
      const d = this.pointer.downAt;
      this.pointer.downAt = null;
      if (d && e.target === canvas && Math.hypot(e.clientX - d.x, e.clientY - d.y) < 5) {
        Object.assign(this.pointer, { x: e.clientX, y: e.clientY });
        this.onAtomClick?.(this._pick());
      }
    });

    // --- measurement overlay: highlighted atoms, dashed path, distance labels ---
    this.selection = [];
    this.onAtomClick = null; // (atomIndex | -1) => void
    this.selGroup = new THREE.Group();
    this.selGroup.renderOrder = 2;
    this.group.add(this.selGroup);
    this.selMat = new THREE.MeshBasicMaterial({ color: 0xffd54a, transparent: true, opacity: 0.38, depthWrite: false });
    this.pathMat = new THREE.LineDashedMaterial({
      color: 0xffd54a, dashSize: 0.18, gapSize: 0.12, depthTest: false, transparent: true,
    });
    this.labelLayer = document.createElement('div');
    this.labelLayer.className = 'measure-labels';
    container.appendChild(this.labelLayer);

    new ResizeObserver(() => this._resize()).observe(container);
    this._resize();
    renderer.setAnimationLoop(() => this._tick());
  }

  // ----------------------------------------------------------------------
  // Public API
  // ----------------------------------------------------------------------

  /** Load a new molecule: re-centres, re-orients and fits the camera. */
  setMolecule(atoms, bonds) {
    this.atoms = atoms;
    this.bonds = bonds;
    this._orient(atoms);
    this._rebuild();
    this.fitView();
  }

  /** Swap coordinates (e.g. next trajectory frame) without moving the camera. */
  setFrame(atoms, bonds) {
    this.atoms = atoms;
    this.bonds = bonds;
    this._rebuild();
  }

  /**
   * Highlight atoms picked for measuring (in click order), draw a dashed path through
   * them and label each segment with its length.
   */
  setSelection(indices) {
    this.selection = indices.slice();
    this._drawSelection();
  }

  _drawSelection() {
    for (const obj of [...this.selGroup.children]) {
      this.selGroup.remove(obj);
      if (obj.isLine) obj.geometry.dispose();
    }
    this.labelLayer.replaceChildren();
    this.segments = [];
    const sel = this.selection.filter((i) => i < this.atoms.length);
    if (!sel.length) return;

    const pts = sel.map((i) => new THREE.Vector3(this.atoms[i].x, this.atoms[i].y, this.atoms[i].z));
    sel.forEach((i, k) => {
      const m = new THREE.Mesh(this._sphereGeo(0), this.selMat);
      m.position.copy(pts[k]);
      m.scale.setScalar((this.radii[i] || 0.3) + 0.3);
      m.renderOrder = 2;
      this.selGroup.add(m);
    });
    if (pts.length > 1) {
      const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), this.pathMat);
      line.computeLineDistances();
      line.renderOrder = 3;
      this.selGroup.add(line);
      for (let k = 1; k < pts.length; k++) {
        const el = document.createElement('div');
        el.className = 'measure-label';
        el.textContent = `${pts[k - 1].distanceTo(pts[k]).toFixed(3)} Å`;
        this.labelLayer.appendChild(el);
        this.segments.push({ mid: pts[k - 1].clone().add(pts[k]).multiplyScalar(0.5), el });
      }
    }
  }

  _placeLabels() {
    if (!this.segments?.length) return;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const v = new THREE.Vector3();
    for (const { mid, el } of this.segments) {
      v.copy(mid);
      this.group.localToWorld(v);
      v.project(this.camera);
      const visible = v.z < 1;
      el.style.display = visible ? '' : 'none';
      if (visible) {
        el.style.left = `${((v.x + 1) / 2) * rect.width}px`;
        el.style.top = `${((1 - v.y) / 2) * rect.height}px`;
      }
    }
  }

  setStyle(style) {
    if (!STYLES[style] || style === this.style) return;
    this.style = style;
    this._rebuild();
    // Bigger atoms (e.g. space-filling) may no longer fit: back the camera off, keeping the rotation.
    const need = this._fitDistance();
    const cam = this.camera;
    const target = this.controls.target;
    if (need > cam.position.distanceTo(target)) {
      cam.position.sub(target).setLength(need).add(target);
      this._setClipPlanes(need);
    }
  }

  /** Excited-state glow sources: [{ x, y, z, weight, color }] in molecule coordinates ([] = off). */
  setGlow(sources) {
    this.glow.setSources(sources);
  }

  /** { sigma (Å), gain } for the glow. */
  setGlowParams(params) {
    this.glow.setParams(params);
  }

  setShowHydrogens(show) {
    this.showHydrogens = show;
    this._rebuild();
  }

  /** Point the camera straight down the molecule's flattest axis and zoom to fit. */
  fitView() {
    const dist = this._fitDistance();
    const cam = this.camera;
    const c = this.controls;
    c.target.set(0, 0, 0);
    cam.position.set(0, 0, dist);
    cam.up.set(0, 1, 0);
    cam.lookAt(0, 0, 0);
    this._setClipPlanes(dist);
    // kill any leftover spin / zoom momentum
    c._movePrev.copy(c._moveCurr);
    c._lastAngle = 0;
    c._zoomStart.copy(c._zoomEnd);
    c._panStart.copy(c._panEnd);
  }

  /** Camera distance at which the whole (visible) molecule fits on screen. */
  _fitDistance() {
    let R = 0;
    for (const i of this.visible) {
      const a = this.atoms[i];
      const d = Math.hypot(a.x - this.center.x, a.y - this.center.y, a.z - this.center.z);
      R = Math.max(R, d + this.radii[i]);
    }
    R = Math.max(R, 1);
    const vfov = THREE.MathUtils.degToRad(this.camera.fov);
    const hfov = 2 * Math.atan(Math.tan(vfov / 2) * this._visibleAspect());
    return (1.08 * R) / Math.sin(Math.min(vfov, hfov) / 2);
  }

  _setClipPlanes(dist) {
    const cam = this.camera;
    cam.near = Math.max(0.01, dist / 200);
    cam.far = dist * 50;
    cam.updateProjectionMatrix();
    this.controls.minDistance = 0.5;
    this.controls.maxDistance = dist * 8;
  }

  // ----------------------------------------------------------------------
  // Scene building
  // ----------------------------------------------------------------------

  _sphereGeo(atomCount) {
    // Fewer triangles per sphere for very large structures.
    const [w, h] = atomCount > 20000 ? [12, 8] : atomCount > 3000 ? [20, 14] : [32, 20];
    const key = `${w}x${h}`;
    return (this.sphereGeos[key] ??= new THREE.SphereGeometry(1, w, h));
  }

  /**
   * Centre the molecule on its centroid and rotate it so its longest extent runs
   * left–right, the next along the screen's vertical, and the thinnest points at
   * the viewer. Flat molecules (benzene, porphyrins...) therefore open face-on.
   */
  _orient(atoms) {
    const n = atoms.length;
    const c = new THREE.Vector3();
    for (const a of atoms) c.add(new THREE.Vector3(a.x, a.y, a.z));
    c.divideScalar(Math.max(n, 1));
    this.center.copy(c);
    this.group.position.copy(c).negate();

    // covariance matrix
    const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const a of atoms) {
      const d = [a.x - c.x, a.y - c.y, a.z - c.z];
      for (let r = 0; r < 3; r++) for (let s = 0; s < 3; s++) C[r][s] += d[r] * d[s];
    }
    const { vectors } = jacobiEigen(C); // sorted by descending eigenvalue
    const e1 = new THREE.Vector3(...vectors[0]);
    const e2 = new THREE.Vector3(...vectors[1]);
    const e3 = new THREE.Vector3().crossVectors(e1, e2); // keep it right-handed
    // Rotation that maps e1 -> x, e2 -> y, e3 -> z is the transpose of the basis matrix.
    const m = new THREE.Matrix4().makeBasis(e1, e2, e3).transpose();
    this.pivot.quaternion.setFromRotationMatrix(m);
  }

  _clear() {
    for (const obj of this.meshes) {
      this.group.remove(obj);
      if (obj.isInstancedMesh) obj.dispose();
      else obj.geometry.dispose(); // line geometry is per-build
    }
    this.meshes = [];
    this.atomMesh = null;
  }

  _rebuild() {
    this._clear();
    const style = STYLES[this.style];
    const atoms = this.atoms;
    const n = atoms.length;

    const show = new Uint8Array(n);
    const visible = [];
    for (let i = 0; i < n; i++) {
      if (this.showHydrogens || atoms[i].symbol !== 'H') {
        show[i] = 1;
        visible.push(i);
      }
    }
    this.visible = visible;

    const els = atoms.map((a) => getElement(a.symbol));
    const radii = (this.radii = new Float32Array(n));
    for (let i = 0; i < n; i++) radii[i] = style.atomRadius(els[i]);

    const m = new THREE.Matrix4();
    const col = new THREE.Color();

    // ---- atoms: one InstancedMesh, one draw call ----
    const atomMesh = new THREE.InstancedMesh(this._sphereGeo(n), this.atomMat, Math.max(1, visible.length));
    atomMesh.count = visible.length;
    visible.forEach((i, k) => {
      const a = atoms[i];
      const r = radii[i];
      m.makeScale(r, r, r).setPosition(a.x, a.y, a.z);
      atomMesh.setMatrixAt(k, m);
      atomMesh.setColorAt(k, col.setHex(els[i].color));
    });
    atomMesh.instanceMatrix.needsUpdate = true;
    if (atomMesh.instanceColor) atomMesh.instanceColor.needsUpdate = true;
    atomMesh.layers.enable(DEPTH_LAYER); // atoms hide the glow behind them
    this._add(atomMesh);
    this.atomMesh = atomMesh;

    // ---- bonds that connect two visible atoms ----
    const vb = [];
    for (let b = 0; b < this.bonds.length; b += 2) {
      const i = this.bonds[b];
      const j = this.bonds[b + 1];
      if (show[i] && show[j]) vb.push(i, j);
    }
    const nb = vb.length / 2;

    if (style.bondRadius > 0 && nb > 0) {
      // Each bond = two half-cylinders, each coloured like the atom it touches.
      const mesh = new THREE.InstancedMesh(this.cylGeo, this.bondMat, nb * 2);
      const A = new THREE.Vector3();
      const B = new THREE.Vector3();
      const mid = new THREE.Vector3();
      const dir = new THREE.Vector3();
      const pos = new THREE.Vector3();
      const scl = new THREE.Vector3();
      const q = new THREE.Quaternion();
      const r = style.bondRadius;
      let k = 0;
      for (let b = 0; b < vb.length; b += 2) {
        const i = vb[b];
        const j = vb[b + 1];
        A.set(atoms[i].x, atoms[i].y, atoms[i].z);
        B.set(atoms[j].x, atoms[j].y, atoms[j].z);
        mid.addVectors(A, B).multiplyScalar(0.5);
        for (const [P, idx] of [[A, i], [B, j]]) {
          dir.subVectors(mid, P);
          const len = dir.length();
          dir.divideScalar(len);
          q.setFromUnitVectors(Y_AXIS, dir);
          pos.addVectors(P, mid).multiplyScalar(0.5);
          scl.set(r, len, r);
          m.compose(pos, q, scl);
          mesh.setMatrixAt(k, m);
          mesh.setColorAt(k, col.setHex(els[idx].color));
          k++;
        }
      }
      mesh.instanceMatrix.needsUpdate = true;
      mesh.instanceColor.needsUpdate = true;
      mesh.layers.enable(DEPTH_LAYER);
      this._add(mesh);
    }

    if (style.lines && nb > 0) {
      const positions = new Float32Array(nb * 12);
      const colors = new Float32Array(nb * 12);
      let p = 0;
      for (let b = 0; b < vb.length; b += 2) {
        const a = atoms[vb[b]];
        const c = atoms[vb[b + 1]];
        const mx = (a.x + c.x) / 2;
        const my = (a.y + c.y) / 2;
        const mz = (a.z + c.z) / 2;
        const ca = col.setHex(els[vb[b]].color).toArray();
        const cc = col.setHex(els[vb[b + 1]].color).toArray();
        positions.set([a.x, a.y, a.z, mx, my, mz, mx, my, mz, c.x, c.y, c.z], p);
        colors.set([...ca, ...ca, ...cc, ...cc], p);
        p += 12;
      }
      const geo = new LineSegmentsGeometry();
      geo.setPositions(positions);
      geo.setColors(colors);
      this._add(new LineSegments2(geo, this.lineMat));
    }

    this._drawSelection(); // atoms may have moved (trajectory) or changed size (style)

    // keep the hover halo on the same atom if it still exists
    if (this.hovered >= n || (this.hovered >= 0 && !show[this.hovered])) this._setHover(-1);
    else if (this.hovered >= 0) this._placeHalo(this.hovered);
    this.pointer.dirty = true;
  }

  _add(obj) {
    this.group.add(obj);
    this.meshes.push(obj);
  }

  // ----------------------------------------------------------------------
  // Hover
  // ----------------------------------------------------------------------

  _pick() {
    if (!this.atomMesh || this.atomMesh.count === 0) return -1;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
      ((this.pointer.x - rect.left) / rect.width) * 2 - 1,
      -((this.pointer.y - rect.top) / rect.height) * 2 + 1,
    );
    this.raycaster.setFromCamera(ndc, this.camera);
    const hit = this.raycaster.intersectObject(this.atomMesh, false)[0];
    return hit ? this.visible[hit.instanceId] : -1;
  }

  _placeHalo(i) {
    const a = this.atoms[i];
    const r = this.radii[i] + 0.25;
    this.halo.position.set(a.x, a.y, a.z);
    this.halo.scale.setScalar(r);
    this.halo.visible = true;
  }

  _setHover(i) {
    this.hovered = i;
    if (i < 0) this.halo.visible = false;
    else this._placeHalo(i);
    this.onHover?.(i, i >= 0 ? this.atoms[i] : null, this.pointer.x, this.pointer.y);
  }

  // ----------------------------------------------------------------------
  // Loop
  // ----------------------------------------------------------------------

  /**
   * Pixels on the left covered by UI (the side panel). The view is shifted so the
   * molecule is centred in the uncovered part of the window instead of behind the panel.
   */
  setInsetLeft(px) {
    this.insetLeft = px;
    this._resize();
  }

  _effectiveInset() {
    const w = this.container.clientWidth;
    // On narrow windows the panel just overlaps; don't squeeze the molecule into a sliver.
    return w - (this.insetLeft ?? 0) >= 320 ? (this.insetLeft ?? 0) : 0;
  }

  _visibleAspect() {
    const w = this.container.clientWidth;
    const h = this.container.clientHeight || 1;
    return (w - this._effectiveInset()) / h;
  }

  _resize() {
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.setViewOffset(w, h, -this._effectiveInset() / 2, 0, w, h);
    this.camera.updateProjectionMatrix();
    this.controls.handleResize();
    this.lineMat.resolution.set(w, h);
    const size = this.renderer.getDrawingBufferSize(new THREE.Vector2());
    this.glow.setSize(size.x, size.y);
  }

  _tick() {
    this.controls.update();
    if (this.pointer.dirty && this.pointer.inside && !this.pointer.down) {
      this.pointer.dirty = false;
      this._setHover(this._pick());
    }
    this.renderer.render(this.scene, this.camera);
    this.glow.render(this.scene, this.camera, Math.min(this.clock.getDelta(), 0.1));
    this._placeLabels();
  }
}

/**
 * Eigen-decomposition of a symmetric 3×3 matrix (cyclic Jacobi).
 * Returns { values, vectors } sorted by descending eigenvalue; vectors[k] is a unit 3-vector.
 */
function jacobiEigen(S) {
  const a = S.map((row) => row.slice());
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    const off = a[0][1] ** 2 + a[0][2] ** 2 + a[1][2] ** 2;
    if (off < 1e-18) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      if (Math.abs(a[p][q]) < 1e-20) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1);
      const s = t * c;
      for (let k = 0; k < 3; k++) {
        const akp = a[k][p];
        const akq = a[k][q];
        a[k][p] = c * akp - s * akq;
        a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = a[p][k];
        const aqk = a[q][k];
        a[p][k] = c * apk - s * aqk;
        a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k][p];
        const vkq = v[k][q];
        v[k][p] = c * vkp - s * vkq;
        v[k][q] = s * vkp + c * vkq;
      }
    }
  }
  const order = [0, 1, 2].sort((i, j) => a[j][j] - a[i][i]);
  return {
    values: order.map((i) => a[i][i]),
    vectors: order.map((i) => [v[0][i], v[1][i], v[2][i]]),
  };
}
