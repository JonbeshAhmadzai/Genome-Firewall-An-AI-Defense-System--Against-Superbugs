const esc = value => String(value ?? '').replace(/[&<>"']/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
const percent = value => value === null || value === undefined || Number.isNaN(Number(value)) ? '—' : `${(Number(value) * 100).toFixed(1)}%`;
const card = (label, value) => `<div class="card"><div class="card-label">${esc(label)}</div><div class="card-value">${esc(value)}</div></div>`;
function comparisonTable(rows) {
  const cols = [['antibiotic','Antibiotic'],['model','Model'],['balanced_accuracy','Balanced accuracy'],['resistant_recall','Resistant recall'],['susceptible_recall','Susceptible recall'],['f1','F1'],['brier','Brier'],['no_call_rate','No-call rate'],['called_accuracy','Called accuracy']];
  if (!rows.length) return '<p class="hint">No benchmark table is available.</p>';
  return `<table><thead><tr>${cols.map(([, label]) => `<th>${label}</th>`).join('')}</tr></thead><tbody>${rows.map(row => `<tr>${cols.map(([key]) => `<td>${key === 'antibiotic' || key === 'model' ? esc(row[key]) : key === 'brier' ? Number(row[key] ?? 0).toFixed(3) : percent(row[key])}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
}
function metricsTable(rows) {
  const cols = [['antibiotic','Antibiotic'],['n_test','Test genomes'],['balanced_accuracy','Balanced accuracy'],['resistant_recall','Resistant recall'],['susceptible_recall','Susceptible recall'],['f1','F1'],['auroc','AUROC'],['pr_auc','PR-AUC'],['brier','Brier'],['no_call_rate','No-call rate']];
  if (!rows.length) return '<p class="hint">No held-out metrics are available.</p>';
  return `<table><thead><tr>${cols.map(([, label]) => `<th>${label}</th>`).join('')}</tr></thead><tbody>${rows.map(row => `<tr>${cols.map(([key]) => `<td>${key === 'antibiotic' ? esc(row[key]) : key === 'n_test' ? Number(row[key] || 0) : key === 'brier' ? Number(row[key] ?? 0).toFixed(3) : percent(row[key])}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
}
async function load() {
  const response = await fetch('/api/overview');
  if (!response.ok) throw new Error('Coverage could not be loaded.');
  const data = await response.json();
  const summary = data.summary || {};
  document.querySelector('#overview-cards').innerHTML = [card('Genomes', Number(summary.genomes || 0).toLocaleString()), card('Drug models', Number(summary.drug_models || 0).toLocaleString()), card('AMR features', Number(summary.amr_features || 0).toLocaleString()), card('Evidence rows', Number(summary.evidence_rows || 0).toLocaleString())].join('');
  const scope = data.scope || {};
  document.querySelector('#scope-status').textContent = scope.status || 'Research scope';
  document.querySelector('#scope-line').innerHTML = `<strong>Species:</strong> ${(scope.species || []).map(esc).join(', ') || 'none'} · <strong>Antibiotics:</strong> ${(scope.antibiotics || []).map(esc).join(', ') || 'none'}`;
  document.querySelector('#model-comparison').innerHTML = comparisonTable(data.model_comparison || []);
  document.querySelector('#held-out-metrics').innerHTML = metricsTable(data.held_out_metrics || []);
  const policy = data.decision_policy || {}, split = data.validation_split || {};
  const steps = (data.pipeline || []).map(step => `<li>${esc(step)}</li>`).join('');
  document.querySelector('#pipeline-info').innerHTML = `<p><strong>Confidence threshold:</strong> ${percent(policy.minimum_confidence)} · <strong>No-call:</strong> ${esc(policy.no_call || '')}</p><p><strong>Validation:</strong> ${percent(split.train)} train · ${percent(split.calibration)} calibration · ${percent(split.test)} grouped test</p><p>${esc(split.method || '')}. ${esc(policy.safety || '')}</p><ol>${steps}</ol>`;
}
load().catch(error => { document.querySelector('#scope-status').textContent = 'Unavailable'; document.querySelector('#scope-line').textContent = error.message; });
