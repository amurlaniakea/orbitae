// Copyright (C) Pedro Sordo Martínez
// SPDX-License-Identifier: AGPL-3.0-or-later
// Prueba funcional del visor: carga web/index.html con los datos de web/data en un navegador
// simulado (jsdom) y un contexto WebGL vacío. Comprueba la lógica y la interfaz, NO el dibujo 3D.
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '../../web');
if (!fs.existsSync(path.join(root, 'data', 'atoms.json')) || !fs.existsSync(path.join(root, 'data', 'molecules.json'))) {
  console.error('Faltan los datos del visor. Genéralos antes con python -m orbitae.export / atoms / hf_compare / molecules (--out web/data).');
  process.exit(2);
}
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const glStub = () => new Proxy({}, { get: (_t, k) => {
  if (k === 'getShaderParameter' || k === 'getProgramParameter') return () => true;
  if (k === 'getShaderInfoLog') return () => '';
  return () => ({});
} });
const ctx2d = new Proxy({}, { get: () => () => {}, set: () => true });
const pageErrors = [];
const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true, url: 'http://localhost/',
  beforeParse(w) {
    w.HTMLCanvasElement.prototype.getContext = function (t) { return t === 'webgl' ? glStub() : ctx2d; };
    w.fetch = async (u) => {
      const p = path.join(root, u);
      if (!fs.existsSync(p)) return { ok: false, status: 404 };
      const b = fs.readFileSync(p);
      return { ok: true, json: async () => JSON.parse(b.toString('utf8')),
               arrayBuffer: async () => b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) };
    };
    w.addEventListener('error', (e) => pageErrors.push(e.message));
  },
});
const w = dom.window, d = w.document, $ = (s) => d.querySelector(s);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const info = () => [...$('#info').children].map((e) => e.textContent).join(' | ');
const num = (s) => s.replace(',', '.');
let fails = 0;
const ok = (c, msg) => { console.log((c ? 'OK     ' : 'FALLA  ') + msg); if (!c) fails++; };
const click = (sel, pred) => {
  const b = [...d.querySelectorAll(sel)].find(pred);
  if (!b) throw new Error('botón no encontrado: ' + sel);
  b.click();
};

(async () => {
  await sleep(500);
  ok($('#msg').style.display !== 'grid', 'arranque sin error visible');
  ok(info().includes('2p') && info().includes('Energía'), 'hidrógeno: información del orbital inicial');
  click('#modes button', (b) => b.dataset.m === 'a'); await sleep(600);
  ok(info().includes('Carbono'), 'átomos: carbono por defecto');
  click('#pt button', (b) => b.title.startsWith('Hierro')); await sleep(300);
  click('#asel button', (b) => b.dataset.k === '3d:-2'); await sleep(300);
  ok(/Subcapa 3d HF/.test(info()) && /Diferencia subcapa \| \+5\d[.,]\d/.test(info()), 'átomos: Fe, subcapa 3d frente a HF (+50 a +59 %)');
  ok(/Diferencia \| \+1[89][.,]\d/.test(info()), 'átomos: Fe, radio medio total (+18 a +19 %, con la configuración forzada)');
  click('#ions button', (b) => b.textContent === 'Fe²⁺'); await sleep(300);
  ok(info().includes('Carga') && info().includes('+2') && info().includes('Fe²⁺'), 'iones: Fe²⁺ con carga +2');
  ok(info().includes('HF frente al neutro'), 'iones: fila frente al átomo neutro');
  const asel = $('#asel').textContent;
  ok(asel.includes('3d') && !asel.includes('4s'), 'iones: Fe²⁺ conserva 3d y pierde 4s');
  click('#modes button', (b) => b.dataset.m === 'm'); await sleep(700);
  ok(info().includes('Agua (H₂O)') && info().includes('HOMO'), 'moléculas: agua con su HOMO');
  ok($('#radsec').hidden === true && $('#mpanel').hidden === false, 'moléculas: panel propio y gráfica radial oculta');
  click('#mols button', (b) => b.textContent === 'CO'); await sleep(600);
  ok(info().includes('Monóxido') && $('#mnote').textContent.includes('carbono'), 'moléculas: CO con la nota de su HOMO');
  click('#msel button', (b) => b.textContent === 'LUMO (1/2)'); await sleep(400);
  ok($('#mnote').textContent.includes('virtuales'), 'moléculas: aviso de orbital virtual');
  ok(d.querySelectorAll('#msel button').length === 7, 'moléculas: 7 orbitales en el CO');
  ok($('#txt').textContent.includes('molécula'), 'moléculas: texto de pestaña propio');
  click('#modes button', (b) => b.dataset.m === 'h'); await sleep(400);
  ok($('#radsec').hidden === false && info().includes('Energía'), 'vuelta al modo hidrógeno');
  ok($('#msg').style.display !== 'grid', 'sin errores visibles al final');
  ok(pageErrors.length === 0, 'sin excepciones en la página' + (pageErrors.length ? ': ' + pageErrors[0] : ''));
  console.log(fails ? `\n${fails} comprobaciones fallaron` : '\nTodas las comprobaciones del visor pasaron');
  w.close();
  process.exit(fails ? 1 : 0);
})().catch((e) => { console.error('EXCEPCIÓN:', e.message); process.exit(2); });
