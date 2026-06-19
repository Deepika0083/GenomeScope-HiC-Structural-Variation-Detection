/* ===============================
   GenomeScope – Clean Version
   =============================== */

// ================= GLOBAL STATE =================
let currentSVs = [];

// ================= FILE LOAD =================
function loadFile(input) {
  if (!input.files[0]) return;
  showToast("File loaded: " + input.files[0].name + " — click Run Analysis");
}

// ================= RUN ANALYSIS =================
async function runAnalysis() {
  const fileInput = document.querySelector("input[type='file']");
  const file = fileInput.files[0];

  if (!file) {
    alert("Please upload a CSV file first!");
    return;
  }

  const formData = new FormData();
  formData.append("file", file);

  // Show loading
  document.getElementById("progress-card").style.display = "block";
  document.getElementById("results-area").style.display = "none";

  try {
    const res = await fetch("http://127.0.0.1:5000/analyze", {
      method: "POST",
      body: formData
    });

    const data = await res.json();
    console.log("Backend result:", data);

    // Hide loading
    document.getElementById("progress-card").style.display = "none";

    // Show results
    showResults(data);

  } catch (err) {
    alert("Backend not connected!");
    console.error(err);
  }
}

// ================= SHOW RESULTS =================
function showResults(data) {
  currentSVs = data;

  document.getElementById("results-area").style.display = "block";
  document.getElementById("metrics-row").style.display = "grid";

  // SV count
  document.getElementById("m-svcount").textContent = data.length;

  // Average confidence
  const avgConf = (
    data.reduce((a, b) => a + b.Confidence, 0) / data.length * 100
  ).toFixed(0) + "%";

  document.getElementById("m-conf").textContent = avgConf;

  // Max size
  const sizes = data
    .filter(s => s.Size_kb)
    .map(s => s.Size_kb / 1000);

  document.getElementById("m-size").textContent =
    sizes.length ? Math.max(...sizes).toFixed(1) + " Mb" : "-";

  // Draw table
  drawTable(data);
}

// ================= DRAW TABLE =================
const typeClass = {
  DEL: "type-del",
  DUP: "type-dup",
  INV: "type-inv",
  TRA: "type-tra"
};

function drawTable(svs) {
  const tbody = document.getElementById("sv-tbody");
  tbody.innerHTML = "";

  svs.forEach(sv => {
    const conf = sv.Confidence || 0;

    const barColor =
      conf > 0.8 ? "#22c55e" :
      conf > 0.6 ? "#f59e0b" :
      "#ef4444";

    const tr = document.createElement("tr");

    tr.innerHTML = `
      <td>${sv.SV_ID || "-"}</td>
      <td><span class="${typeClass[sv.Type] || "type-del"}">${sv.Type}</span></td>
      <td>${sv.Chromosome || "-"}</td>
      <td>${sv.Start || "-"}</td>
      <td>${sv.End || "-"}</td>
      <td>${sv.Size_kb ? (sv.Size_kb / 1000).toFixed(2) + " Mb" : "-"}</td>
      <td>
        <div style="display:flex;align-items:center;gap:6px;">
          <div style="width:80px;height:6px;background:#222;border-radius:4px;">
            <div style="width:${(conf * 100).toFixed(0)}%;height:100%;background:${barColor};"></div>
          </div>
          <span>${(conf * 100).toFixed(0)}%</span>
        </div>
      </td>
      <td>${sv.Disease || "-"}</td>
    `;

    tbody.appendChild(tr);
  });
}

// ================= EXPORT CSV =================
function exportReport() {
  if (!currentSVs.length) {
    showToast("Run analysis first!");
    return;
  }

  let csv = "SV_ID,Type,Chromosome,Start,End,Size_kb,Confidence,Disease\n";

  currentSVs.forEach(sv => {
    csv += `${sv.SV_ID},${sv.Type},${sv.Chromosome},${sv.Start},${sv.End},${sv.Size_kb},${sv.Confidence},${sv.Disease}\n`;
  });

  const blob = new Blob([csv], { type: "text/csv" });
  const a = document.createElement("a");

  a.href = URL.createObjectURL(blob);
  a.download = "genomescope_report.csv";
  a.click();
}

// ================= TOAST =================
function showToast(msg) {
  const t = document.createElement("div");

  t.style.cssText = `
    position:fixed;
    bottom:20px;
    right:20px;
    background:#111;
    color:#fff;
    padding:10px 15px;
    border-radius:8px;
    z-index:9999;
  `;

  t.textContent = msg;
  document.body.appendChild(t);

  setTimeout(() => t.remove(), 3000);
}

// ================= DRAG & DROP =================
const dropZone = document.getElementById("drop-zone");

if (dropZone) {
  dropZone.addEventListener("dragover", e => {
    e.preventDefault();
    dropZone.classList.add("drag");
  });

  dropZone.addEventListener("dragleave", () => {
    dropZone.classList.remove("drag");
  });

  dropZone.addEventListener("drop", e => {
    e.preventDefault();
    dropZone.classList.remove("drag");

    const file = e.dataTransfer.files[0];
    if (file) {
      document.querySelector("input[type='file']").files = e.dataTransfer.files;
      showToast("File loaded: " + file.name);
    }
  });
}