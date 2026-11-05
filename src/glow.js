// Excited-state glow: a soft, emissive Gaussian cloud drawn over part of a molecule.
//
// Each source (usually an atom of the excited chromophore) is a 3D Gaussian of light. Drawn as a camera-facing
// sprite with a 2D Gaussian falloff and added up, that is exactly the line-of-sight integral of the 3D cloud, so
// overlapping sources merge into one smooth glow that follows the molecule's shape.
//
// Passes (only while something glows):
//   1. glow target: depth of atoms + bonds (so atoms in front hide the glow), then the sprites, added up in HDR;
//   2. a wide blur of that at quarter resolution (bloom on the glow only, the rest of the molecule stays crisp);
//   3. composite onto the screen with a saturating curve, 1 - exp(-I): faint parts keep the glow colour, dense
//      parts run toward white, like a hot core.

import * as THREE from 'three';
import { FullScreenQuad } from 'three/addons/postprocessing/Pass.js';

export const GLOW_LAYER = 1; // sprites
export const DEPTH_LAYER = 2; // atoms + bonds, for occluding the glow

const spriteVert = /* glsl */ `
  attribute vec3 iCenter;
  attribute vec3 iColor;
  attribute float iWeight;
  attribute float iPhase;
  uniform float uSigma;
  uniform float uTime;
  uniform float uPush;
  uniform float uGain;
  varying vec2 vR;
  varying vec3 vColor;
  void main() {
    vec4 mv = modelViewMatrix * vec4(iCenter, 1.0);
    mv.z += uPush;                       // nudge toward the camera
    mv.xy += position.xy * 3.0 * uSigma; // sprite spans +-3 sigma
    vR = position.xy * 3.0;              // distance from the centre in units of sigma
    // breathing pulse + a little per-atom shimmer
    float pulse = 1.0 + 0.16 * sin(6.2831853 * uTime / 1.9) + 0.07 * sin(4.7 * uTime + iPhase);
    vColor = iColor * iWeight * uGain * pulse;
    gl_Position = projectionMatrix * mv;
  }
`;

const spriteFrag = /* glsl */ `
  varying vec2 vR;
  varying vec3 vColor;
  void main() {
    float g = exp(-0.5 * dot(vR, vR));
    gl_FragColor = vec4(vColor * g, 1.0);
  }
`;

const quadVert = /* glsl */ `
  varying vec2 vUv;
  void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }
`;

const blurFrag = /* glsl */ `
  uniform sampler2D tSrc;
  uniform vec2 uStep;
  varying vec2 vUv;
  void main() {
    vec3 c = texture2D(tSrc, vUv).rgb * 0.2270270270;
    c += (texture2D(tSrc, vUv + uStep * 1.3846153846).rgb + texture2D(tSrc, vUv - uStep * 1.3846153846).rgb) * 0.3162162162;
    c += (texture2D(tSrc, vUv + uStep * 3.2307692308).rgb + texture2D(tSrc, vUv - uStep * 3.2307692308).rgb) * 0.0702702703;
    gl_FragColor = vec4(c, 1.0);
  }
`;

const compositeFrag = /* glsl */ `
  uniform sampler2D tGlow;
  uniform sampler2D tBloom;
  uniform float uBloom;
  varying vec2 vUv;
  void main() {
    vec3 I = texture2D(tGlow, vUv).rgb + uBloom * texture2D(tBloom, vUv).rgb;
    vec3 c = 1.0 - exp(-I);
    gl_FragColor = vec4(c, max(c.r, max(c.g, c.b)));
  }
`;

export class GlowLayer {
  constructor(renderer) {
    this.renderer = renderer;
    this.count = 0;
    this.time = 0;

    const quad = new THREE.PlaneGeometry(2, 2);
    this.geo = new THREE.InstancedBufferGeometry();
    this.geo.index = quad.index;
    this.geo.setAttribute('position', quad.getAttribute('position'));
    this.spriteMat = new THREE.ShaderMaterial({
      uniforms: {
        uSigma: { value: 1.2 },
        uTime: { value: 0 },
        uPush: { value: 0.15 }, // Å toward the camera: mostly behind the atoms (backlit), a slight tint over them
        uGain: { value: 1 },
      },
      vertexShader: spriteVert,
      fragmentShader: spriteFrag,
      blending: THREE.AdditiveBlending,
      transparent: true,
      depthWrite: false,
      depthTest: true,
    });
    this.mesh = new THREE.Mesh(this.geo, this.spriteMat);
    this.mesh.frustumCulled = false;
    this.mesh.layers.set(GLOW_LAYER);

    this.depthMat = new THREE.MeshBasicMaterial({ colorWrite: false });

    const rtOpts = { type: THREE.HalfFloatType, depthBuffer: false };
    this.rt = new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType, depthBuffer: true });
    this.rtA = new THREE.WebGLRenderTarget(1, 1, rtOpts);
    this.rtB = new THREE.WebGLRenderTarget(1, 1, rtOpts);

    this.blurMat = new THREE.ShaderMaterial({
      uniforms: { tSrc: { value: null }, uStep: { value: new THREE.Vector2() } },
      vertexShader: quadVert,
      fragmentShader: blurFrag,
      depthTest: false,
      depthWrite: false,
    });
    this.compositeMat = new THREE.ShaderMaterial({
      uniforms: { tGlow: { value: this.rt.texture }, tBloom: { value: this.rtB.texture }, uBloom: { value: 1.6 } },
      vertexShader: quadVert,
      fragmentShader: compositeFrag,
      depthTest: false,
      depthWrite: false,
      transparent: true,
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneMinusSrcColorFactor, // "screen": brightens, never over-exposes what's underneath
      blendSrcAlpha: THREE.OneFactor,
      blendDstAlpha: THREE.OneMinusSrcAlphaFactor,
    });
    this.fsq = new FullScreenQuad(null);
    this._clear = new THREE.Color();
  }

  /** sources: [{ x, y, z, weight, color: THREE.Color }] in molecule coordinates. */
  setSources(sources) {
    const n = sources.length;
    const center = new Float32Array(n * 3);
    const color = new Float32Array(n * 3);
    const weight = new Float32Array(n);
    const phase = new Float32Array(n);
    sources.forEach((s, i) => {
      center.set([s.x, s.y, s.z], i * 3);
      color.set([s.color.r, s.color.g, s.color.b], i * 3);
      weight[i] = s.weight;
      phase[i] = (i * 2.399963) % 6.2831853; // golden-angle spread: neighbours shimmer out of step
    });
    for (const name of ['iCenter', 'iColor', 'iWeight', 'iPhase']) this.geo.deleteAttribute(name);
    this.geo.setAttribute('iCenter', new THREE.InstancedBufferAttribute(center, 3));
    this.geo.setAttribute('iColor', new THREE.InstancedBufferAttribute(color, 3));
    this.geo.setAttribute('iWeight', new THREE.InstancedBufferAttribute(weight, 1));
    this.geo.setAttribute('iPhase', new THREE.InstancedBufferAttribute(phase, 1));
    this.geo.instanceCount = n;
    this.count = n;
  }

  /** sigma: Gaussian width (Å); gain: brightness multiplier. */
  setParams({ sigma, gain }) {
    if (sigma !== undefined) this.spriteMat.uniforms.uSigma.value = sigma;
    if (gain !== undefined) this.spriteMat.uniforms.uGain.value = gain;
  }

  setSize(w, h) {
    this.rt.setSize(w, h);
    const qw = Math.max(1, Math.round(w / 4));
    const qh = Math.max(1, Math.round(h / 4));
    this.rtA.setSize(qw, qh);
    this.rtB.setSize(qw, qh);
  }

  /** Draw the glow over what is already on screen. Call after the normal scene render. */
  render(scene, camera, dt) {
    if (!this.count) return;
    this.time += dt;
    this.spriteMat.uniforms.uTime.value = this.time;
    const r = this.renderer;
    const autoClear = r.autoClear;
    const clearAlpha = r.getClearAlpha();
    r.getClearColor(this._clear);
    const mask = camera.layers.mask;

    // 1. depth of atoms + bonds, then the glow sprites on top
    r.setClearColor(0x000000, 0);
    r.setRenderTarget(this.rt);
    r.clear(true, true, false);
    r.autoClear = false;
    camera.layers.set(DEPTH_LAYER);
    scene.overrideMaterial = this.depthMat;
    r.render(scene, camera);
    scene.overrideMaterial = null;
    camera.layers.set(GLOW_LAYER);
    r.render(scene, camera);
    camera.layers.mask = mask;

    // 2. bloom: separable blur at quarter resolution
    this.fsq.material = this.blurMat;
    const u = this.blurMat.uniforms;
    u.tSrc.value = this.rt.texture;
    u.uStep.value.set(1.6 / this.rtA.width, 0);
    r.setRenderTarget(this.rtA);
    this.fsq.render(r);
    u.tSrc.value = this.rtA.texture;
    u.uStep.value.set(0, 1.6 / this.rtA.height);
    r.setRenderTarget(this.rtB);
    this.fsq.render(r);

    // 3. composite onto the screen
    r.setRenderTarget(null);
    this.fsq.material = this.compositeMat;
    this.fsq.render(r);

    r.autoClear = autoClear;
    r.setClearColor(this._clear, clearAlpha);
  }
}
