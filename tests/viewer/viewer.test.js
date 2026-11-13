// Copyright (C) Pedro Sordo Martínez
// SPDX-License-Identifier: AGPL-3.0-or-later
// Prueba funcional del visor: carga web/index.html con los datos de web/data en un navegador
// simulado (jsdom) y un contexto WebGL vacío. Comprueba la lógica y la interfaz (en español e
// inglés), NO el dibujo 3D.
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '../../web');
for (const f of ['atoms.json', 'molecules.json', 'index.json', 'hf_compare.json']) {
  if (!fs.existsSync(path.join(root, 'data', f))) {
    console.error('Faltan los datos del visor (' + f + '). Genéralos con python -m orbitae.export / atoms / hf_compare / molecules (--out web/data).');
    process.exit(2);
  }
}
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const glStub = () => new Proxy({}, { get: (_t, k) => {
  if (k === 'getShaderParameter' || k === 'getProgramParameter') return () => true;
  if (k === 'getShaderInfoLog') return () => '';
  return () => ({});
} });
const ctx2d = new Proxy({}, { get: () => () => {}, set: () => true });

function boot(navLang) {
  const pageErrors = [];
  const dom = new JSDOM(html, {
    runScripts: 'dangerously', pretendToBeVisual: true, url: 'http://localhost/',
    beforeParse(w) {
      Object.defineProperty(w.navigator, 'language', { value: navLang, configurable: true });
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
  return { w: dom.window, pageErrors };
}
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
let fails = 0;
const ok = (c, msg) => { console.log((c ? 'OK     ' : 'FALLA  ') + msg); if (!c) fails++; };

(async () => {
  // --- idioma por defecto según el navegador
  const en0 = boot('en-US'); await sleep(500);
  const d0 = en0.w.document;
  ok(d0.querySelector('.sub').textContent.includes('Atomic orbitals'), 'navegador en inglés: la interfaz arranca en inglés');
  ok(d0.documentElement.lang === 'en' && d0.title.includes('Orbital viewer'), 'navegador en inglés: título e idioma del documento');
  en0.w.close();

  // --- sesión principal en español
  const { w, pageErrors } = boot('es-ES'); await sleep(500);
  const d = w.document, $ = (s) => d.querySelector(s);
  const info = () => [...$('#info').children].map((e) => e.textContent).join(' | ');
  const click = (sel, pred) => {
    const b = [...d.querySelectorAll(sel)].find(pred);
    if (!b) throw new Error('botón no encontrado: ' + sel);
    b.click();
  };
  ok($('#msg').style.display !== 'grid', 'arranque sin error visible');
  ok($('.sub').textContent.includes('rigor verificable'), 'navegador en español: la interfaz arranca en español');
  ok(info().includes('2p') && info().includes('Energía'), 'hidrógeno: información del orbital inicial');
  click('#modes button', (b) => b.dataset.m === 'a'); await sleep(600);
  ok(info().includes('Carbono'), 'átomos: carbono por defecto');
  click('#pt button', (b) => b.title.startsWith('Hierro')); await sleep(300);
  click('#asel button', (b) => b.dataset.k === '3d:-2'); await sleep(300);
  ok(/Subcapa 3d HF/.test(info()) && /Diferencia subcapa \| \+5\d[.,]\d/.test(info()), 'átomos: Fe, subcapa 3d frente a HF (+50 a +59 %)');
  ok(/Diferencia \| \+1[89][.,]\d/.test(info()), 'átomos: Fe, radio medio total (+18 a +19 %, con la configuración forzada)');
  ok(/\+56,\d/.test(info()), 'español: decimales con coma');
  click('#ions button', (b) => b.textContent === 'Fe²⁺'); await sleep(300);
  ok(info().includes('Carga') && info().includes('+2') && info().includes('Fe²⁺'), 'iones: Fe²⁺ con carga +2');
  ok(info().includes('HF frente al neutro'), 'iones: fila frente al átomo neutro');
  const asel = $('#asel').textContent;
  ok(asel.includes('3d') && !asel.includes('4s'), 'iones: Fe²⁺ conserva 3d y pierde 4s');

  // --- cambio a inglés estando en un orbital concreto del hierro
  click('#pt button', (b) => b.title.startsWith('Hierro')); await sleep(300);
  click('#asel button', (b) => b.dataset.k === '3d:-2'); await sleep(300);
  click('#langs button', (b) => b.dataset.l === 'en'); await sleep(900);
  ok($('.sub').textContent.includes('Atomic orbitals') && d.title.includes('Orbital viewer'), 'inglés: cabecera y título');
  ok(info().includes('Element') && info().includes('Iron (Fe)') && info().includes('Effective Z'), 'inglés: etiquetas y nombre del elemento -> ' + info().slice(0, 70));
  ok(/Subshell 3d, HF/.test(info()) && /Subshell difference \| \+56\.\d/.test(info()), 'inglés: subcapa 3d frente a HF y decimales con punto');
  ok($('#bio').textContent.startsWith('Biomedical relevance: Center of hemoglobin'), 'inglés: relevancia biomédica del hierro');
  ok(d.querySelector('#pt button[title^="Iron"]') !== null, 'inglés: la tabla periódica se reconstruye en inglés');
  ok([...d.querySelectorAll('#modes button')].map((b) => b.textContent).join('|') === 'Exact hydrogen|Atoms and ions (Slater)|Molecules (HF)', 'inglés: botones de modo');
  ok($('#txt').textContent.includes('Each atom distributes its electrons'), 'inglés: texto de la pestaña de estudiante');
  ok($('#asel').textContent.includes('Total cloud') && d.querySelector('#asel button[aria-pressed="true"]').dataset.k === '3d:-2', 'inglés: se conserva el orbital elegido');
  ok(w.localStorage.getItem('orbitae-lang') === 'en', 'inglés: la elección se recuerda');
  click('#modes button', (b) => b.dataset.m === 'm'); await sleep(700);
  ok(info().includes('Water (H₂O)') && info().includes('Orbital energy'), 'inglés, moléculas: agua y etiquetas');
  ok($('#mnote').textContent.includes('Lone pair of oxygen'), 'inglés, moléculas: nota del HOMO del agua');
  click('#mols button', (b) => b.textContent === 'CO'); await sleep(600);
  click('#msel button', (b) => b.textContent === 'LUMO (1/2)'); await sleep(400);
  ok($('#mnote').textContent.includes('back-donation') && $('#mnote').textContent.includes('HF virtual orbitals'), 'inglés, moléculas: nota del LUMO y aviso de orbital virtual');
  ok($('#bio').textContent.includes('Toxic gas'), 'inglés, moléculas: relevancia biomédica del CO');
  click('#langs button', (b) => b.dataset.l === 'es'); await sleep(900);
  ok(info().includes('Monóxido') && $('#mnote').textContent.includes('retrodonación') && $('#mnote').textContent.includes('virtuales'), 'vuelta a español en moléculas, con el LUMO conservado');
  ok(d.querySelectorAll('#msel button').length === 7, 'moléculas: 7 orbitales en el CO');
  ok($('#radsec').hidden === true && $('#mpanel').hidden === false, 'moléculas: panel propio y gráfica radial oculta');
  click('#modes button', (b) => b.dataset.m === 'h'); await sleep(400);
  ok($('#radsec').hidden === false && info().includes('Energía'), 'vuelta al modo hidrógeno en español');
  ok($('#msg').style.display !== 'grid', 'sin errores visibles al final');
  ok(pageErrors.length === 0, 'sin excepciones en la página' + (pageErrors.length ? ': ' + pageErrors[0] : ''));
  console.log(fails ? `\n${fails} comprobaciones fallaron` : '\nTodas las comprobaciones del visor pasaron');
  w.close();
  process.exit(fails ? 1 : 0);
})().catch((e) => { console.error('EXCEPCIÓN:', e.message); process.exit(2); });
