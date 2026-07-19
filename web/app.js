const state = { response: null };
const speciesEl = document.querySelector('#species');
const fileEl = document.querySelector('#fasta');
const explainEl = document.querySelector('#explain');
const runEl = document.querySelector('#run');
const messageEl = document.querySelector('#message');
const resultsEl = document.querySelector('#results');

function showMessage(text, error = false) {
  messageEl.textContent = text;
  messageEl.classList.toggle('hidden', !text);
  messageEl.classList.toggle('error', error);
  messageEl.style.background = error ? '#ffe8ea' : '#fff2e2';
  messageEl.style.color = error ? '#932f3b' : '#8b5315';
}

function card(label, value) {
  return `<div class="card"><div class="card-label">${esc(label)}</div><div class="card-value">${esc(value)}</div></div>`;
}

function esc(value) {
  return String(value ?? '').replace(/[&<>"']/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
}

function percent(value) {
  return value === null || value === undefined || Number.isNaN(Number(value)) ? '—' : `${(Number(value) * 100).toFixed(1)}%`;
}

function renderMarkdown(source) {
  const safe = esc(source || '');
  const inline = value => value.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>').replace(/\*(.+?)\*/g, '<em>$1</em>');
  const lines = safe.split(/\r?\n/);
  const output = [];
  let listOpen = false;
  const closeList = () => { if (listOpen) { output.push('</ul>'); listOpen = false; } };
  lines.forEach(line => {
    if (/^###\s+/.test(line)) { closeList(); output.push(`<h3>${inline(line.replace(/^###\s+/, ''))}</h3>`); return; }
    if (/^####\s+/.test(line)) { closeList(); output.push(`<h4>${inline(line.replace(/^####\s+/, ''))}</h4>`); return; }
    if (/^\s*[-*]\s+/.test(line)) {
      if (!listOpen) { output.push('<ul>'); listOpen = true; }
      output.push(`<li>${inline(line.replace(/^\s*[-*]\s+/, ''))}</li>`); return;
    }
    closeList();
    if (/^\s*\d+\.\s+/.test(line)) {
      const number = line.match(/^\s*(\d+)\./)[1];
      const copy = line.replace(/^\s*\d+\.\s+/, '');
      output.push(`<p class="numbered"><strong>${number}.</strong> ${inline(copy)}</p>`);
      return;
    }
    if (line.trim()) output.push(`<p>${inline(line)}</p>`);
  });
  closeList();
  return output.join('');
}

function comparisonTable(rows) {
  if (!rows || !rows.length) return '<p class="hint">No benchmark table is available.</p>';
  const columns = [
    ['antibiotic', 'Antibiotic'], ['model', 'Model'], ['balanced_accuracy', 'Balanced accuracy'],
    ['resistant_recall', 'Resistant recall'], ['susceptible_recall', 'Susceptible recall'],
    ['f1', 'F1'], ['brier', 'Brier'], ['no_call_rate', 'No-call rate'], ['called_accuracy', 'Called accuracy'],
  ];
  return `<table><thead><tr>${columns.map(([, label]) => `<th>${label}</th>`).join('')}</tr></thead><tbody>` +
    rows.map(row => `<tr>${columns.map(([key]) => `<td>${key === 'antibiotic' || key === 'model' ? esc(row[key]) : (key === 'brier' ? (row[key] == null ? '—' : Number(row[key]).toFixed(3)) : percent(row[key]))}</td>`).join('')}</tr>`).join('') +
    '</tbody></table>';
}

function heldOutTable(rows) {
  if (!rows || !rows.length) return '<p class="hint">No held-out metrics are available.</p>';
  const columns = [
    ['antibiotic', 'Antibiotic'], ['n_test', 'Test genomes'], ['balanced_accuracy', 'Balanced accuracy'],
    ['resistant_recall', 'Resistant recall'], ['susceptible_recall', 'Susceptible recall'],
    ['f1', 'F1'], ['auroc', 'AUROC'], ['pr_auc', 'PR-AUC'], ['brier', 'Brier'], ['no_call_rate', 'No-call rate'],
  ];
  return `<table><thead><tr>${columns.map(([, label]) => `<th>${label}</th>`).join('')}</tr></thead><tbody>` +
    rows.map(row => `<tr>${columns.map(([key]) => `<td>${key === 'antibiotic' ? esc(row[key]) : key === 'n_test' ? Number(row[key] || 0) : key === 'brier' ? (row[key] == null ? '—' : Number(row[key]).toFixed(3)) : percent(row[key])}</td>`).join('')}</tr>`).join('') +
    '</tbody></table>';
}

function renderOverview(data) {
  const summary = data.summary || {};
  document.querySelector('#overview-cards').innerHTML = [
    card('Genomes', Number(summary.genomes || 0).toLocaleString()),
    card('Drug models', Number(summary.drug_models || 0).toLocaleString()),
    card('AMR features', Number(summary.amr_features || 0).toLocaleString()),
    card('Evidence rows', Number(summary.evidence_rows || 0).toLocaleString()),
  ].join('');
  const scope = data.scope || {};
  document.querySelector('#scope-status').textContent = scope.status || 'Research scope';
  document.querySelector('#scope-line').innerHTML = `<strong>Species:</strong> ${(scope.species || []).map(esc).join(', ') || 'none'} · <strong>Antibiotics:</strong> ${(scope.antibiotics || []).map(esc).join(', ') || 'none'}`;
  const comparisonEl = document.querySelector('#model-comparison');
  const heldOutEl = document.querySelector('#held-out-metrics');
  if (comparisonEl) comparisonEl.innerHTML = comparisonTable(data.model_comparison || []);
  if (heldOutEl) heldOutEl.innerHTML = heldOutTable(data.held_out_metrics || []);
  const policy = data.decision_policy || {};
  const split = data.validation_split || {};
  const steps = (data.pipeline || []).map(step => `<li>${esc(step)}</li>`).join('');
  const pipelineEl = document.querySelector('#pipeline-info');
  if (pipelineEl) pipelineEl.innerHTML = `<p><strong>Confidence threshold:</strong> ${percent(policy.minimum_confidence)} · <strong>No-call:</strong> ${esc(policy.no_call || '')}</p><p><strong>Validation:</strong> ${percent(split.train)} train · ${percent(split.calibration)} calibration · ${percent(split.test)} grouped test</p><p>${esc(split.method || '')}. ${esc(policy.safety || '')}</p><ol>${steps}</ol>`;
}

function render(response) {
  state.response = response;
  const qc = response.qc || {};
  document.querySelector('#qc-cards').innerHTML = [
    card('Contigs', Number(qc.contigs || 0).toLocaleString()),
    card('Assembly size', `${Number(qc.total_bases || 0).toLocaleString()} bp`),
    card('N50', `${Number(qc.n50 || 0).toLocaleString()} bp`),
    card('Ambiguous bases', `${(Number(qc.ambiguous_fraction || 0) * 100).toFixed(2)}%`),
  ].join('');
  document.querySelector('#decisions').innerHTML = (response.predictions || []).map(row => {
    const decision = row.prediction || 'no-call';
    const modelDecision = row.prediction_before_target_gate;
    const normalized = decision.toLowerCase();
    const klass = normalized.includes('resist') || normalized.includes('fail') ? 'fail' : normalized.includes('suscept') || normalized.includes('work') ? 'work' : '';
    const gateNote = modelDecision && modelDecision !== decision ? `<div class="decision-gate">Safety gate changed model call from <strong>${esc(modelDecision)}</strong>.</div>` : '';
    return `<div class="decision ${klass}"><div class="decision-title">${esc(row.antibiotic)}: ${esc(decision)}</div>` +
      `<div class="decision-meta">Confidence ${(Number(row.confidence || 0) * 100).toFixed(0)}% · target ${esc(row.target_status || 'unknown')}</div>` +
      `<div class="decision-meta">${esc(row.target_gate_reason || '')}</div>${gateNote}</div>`;
  }).join('');
  const evidence = response.amr_evidence || [];
  document.querySelector('#evidence').innerHTML = evidence.length ?
    `<table><thead><tr><th>Feature</th><th>Gene</th><th>Class</th><th>Subtype</th></tr></thead><tbody>` +
    evidence.map(row => `<tr><td>${esc(row.feature)}</td><td>${esc(row.gene_symbol)}</td><td>${esc(row.amr_class)}</td><td>${esc(row.subtype)}</td></tr>`).join('') +
    '</tbody></table>' : '<div class="evidence-empty"><span class="evidence-empty-icon">i</span><div><strong>No reportable AMR determinant detected</strong><span>AMRFinderPlus found no matching resistance gene or mutation in this input. This does not prove susceptibility; interpret it with the model confidence and laboratory AST.</span></div></div>';
  const explanationPanel = document.querySelector('#explanation-panel');
  explanationPanel.classList.toggle('hidden', !response.llm_explanation);
  document.querySelector('#explanation').innerHTML = renderMarkdown(response.llm_explanation || '');
  resultsEl.classList.remove('hidden');
}

async function loadConfig() {
  const [configResponse, overviewResponse] = await Promise.all([fetch('/api/config'), fetch('/api/overview')]);
  if (!configResponse.ok || !overviewResponse.ok) throw new Error('The API is online but the overview could not be loaded.');
  const config = await configResponse.json();
  const overview = await overviewResponse.json();
  (config.species || []).forEach(species => speciesEl.add(new Option(species, species)));
  renderOverview(overview);
  document.querySelector('#api-status').textContent = config.llm_available ? 'API online · GPT explanations ready' : 'API online';
}

runEl.addEventListener('click', async () => {
  if (!fileEl.files.length) return showMessage('Choose one FASTA file first.', true);
  runEl.disabled = true; runEl.textContent = 'Running FASTA QC and models…'; showMessage('');
  const body = new FormData();
  body.append('file', fileEl.files[0]); body.append('species', speciesEl.value); body.append('explain', explainEl.checked ? 'true' : 'false');
  try {
    const response = await fetch('/api/predict', { method: 'POST', body });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || 'Prediction failed');
    render(data);
  } catch (error) { showMessage(error.message, true); }
  finally { runEl.disabled = false; runEl.textContent = 'Run genome analysis'; }
});

fileEl.addEventListener('change', () => {
  const file = fileEl.files[0];
  document.querySelector('#file-name').textContent = file ? `${file.name} · ${(file.size / 1024 / 1024).toFixed(2)} MB` : '.fa, .fna, .fasta, or .gz';
});

const uploadZone = document.querySelector('.upload-zone');
['dragenter', 'dragover'].forEach(eventName => uploadZone.addEventListener(eventName, event => {
  event.preventDefault(); uploadZone.classList.add('dragging');
}));
['dragleave', 'drop'].forEach(eventName => uploadZone.addEventListener(eventName, event => {
  event.preventDefault(); uploadZone.classList.remove('dragging');
}));
uploadZone.addEventListener('drop', event => {
  const files = event.dataTransfer.files;
  if (!files.length) return;
  try { fileEl.files = files; fileEl.dispatchEvent(new Event('change')); } catch (_) { showMessage('Use the file picker to select this genome.', true); }
});

document.querySelector('#download').addEventListener('click', () => {
  if (!state.response) return;
  const blob = new Blob([JSON.stringify(state.response, null, 2)], { type: 'application/json' });
  const link = document.createElement('a'); link.href = URL.createObjectURL(blob); link.download = 'genome_firewall_response.json'; link.click(); URL.revokeObjectURL(link.href);
});

loadConfig().catch(error => { document.querySelector('#api-status').textContent = 'API unavailable'; showMessage(error.message, true); });
