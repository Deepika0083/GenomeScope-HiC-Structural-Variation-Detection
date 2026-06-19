# GenomeScope
### Structural Variation Detection from 3D Genome (Hi-C) Data
**Mini Project — Computational Biology / Genomics**

---

## What this project does

We analyse Hi-C contact matrices — grids that show how often each pair of genomic
regions physically touched inside a cell — and automatically flag regions that look
structurally abnormal.

Those abnormal regions are structural variations (SVs): deletions, duplications,
inversions, or translocations that can cause cancer or genetic disorders.

---

## Project structure

```
genomescope/
│
├── index.html              ← Open this in a browser — the full website
│
├── css/
│   └── style.css           ← All the visual styling
│
├── js/
│   └── app.js              ← All the browser-side logic (heatmaps, tables, etc.)
│
├── python/
│   └── sv_detector.py      ← Python pipeline (runs on real data)
│
├── data/
│   └── sample_matrix_normal.csv   ← A small test matrix you can load
│
└── README.md               ← This file
```

---

## How to run the website (VS Code)

**Option 1 — Just double-click (simplest)**
```
Open index.html directly in Chrome or Edge.
Everything works without any server.
```

**Option 2 — VS Code Live Server (recommended)**
1. Install the "Live Server" extension in VS Code
   (search "Live Server" by Ritwick Dey in Extensions tab)
2. Right-click `index.html` → "Open with Live Server"
3. Browser opens automatically at `http://127.0.0.1:5500`

---

## How to run the Python script

**Install dependencies (one time only)**
```bash
pip install numpy scipy scikit-learn matplotlib pandas
```

**Run on built-in demo data**
```bash
cd python
python sv_detector.py
```
This generates a deletion pattern, detects it, prints a report, and saves:
- `sv_report.csv`   — the detected SVs as a table
- `genomescope_output.png` — a 3-panel heatmap figure

**Run with a real matrix file**
```bash
python sv_detector.py --matrix ../data/sample_matrix_normal.csv
```

**All options**
```
--matrix      path to your CSV contact matrix file
--sv_type     demo SV type if no matrix: del / dup / inv / tra  (default: del)
--output      where to save the CSV report    (default: sv_report.csv)
--plot        where to save the PNG figure    (default: genomescope_output.png)
--sensitivity z-score cutoff for anomaly scan (default: 3.5, lower = more sensitive)
```

**Examples**
```bash
# Test each SV type
python sv_detector.py --sv_type del
python sv_detector.py --sv_type tra
python sv_detector.py --sv_type inv

# Run on your own matrix with lower sensitivity
python sv_detector.py --matrix my_hic.csv --sensitivity 2.5 --output my_svs.csv
```

---

## Input matrix format

The script expects a plain CSV where:
- Each row and column represents one genomic bin (default 10 kb)
- Cell value = contact frequency between those two bins
- No header row, no row labels

Example (5×5 matrix, 5 bins):
```
8.5, 3.2, 1.1, 0.4, 0.2
3.2, 9.1, 3.4, 1.2, 0.5
1.1, 3.4, 8.8, 3.3, 1.0
0.4, 1.2, 3.3, 9.0, 3.1
0.2, 0.5, 1.0, 3.1, 8.6
```

For real Hi-C data, export a region from:
- **Juicebox** (File → Export → contact matrix as text)
- **cooler**: `cooler dump -t pixels contacts.cool > matrix.csv`
- **HiC-Pro**: the `matrix/` output folder contains `.matrix` files

---

## How the detection algorithm works

```
1. Load matrix
2. O/E normalize  →  divide each pixel by the expected frequency at that distance
3. Extract features  →  for each off-diagonal window: mean, max, std, z-score
4. Isolation Forest  →  unsupervised ML model flags statistical outliers
5. Z-score scan  →  independently flags pixels above a threshold
6. Cluster  →  merge nearby flagged pixels into a single SV event
7. Classify  →  rule-based type assignment (DEL/DUP/INV/TRA) from block geometry
8. Report  →  print table + save CSV + save heatmap figure
```

---

## Public datasets used in the project

| Dataset | Cell line | Source |
|---------|-----------|--------|
| K562 Hi-C (leukemia) | K562 | GEO: GSE63525 |
| GM12878 Hi-C (normal) | GM12878 | GEO: GSE63525 |
| T47D Hi-C (breast cancer) | T47D | GEO: GSE176679 |
| Charm simulated SVs | synthetic | Figshare: 10.6084/m9.figshare.28269734 |
| Reference genome | GRCh38 | Ensembl release 107 |

---

## What each file contains

| File | What's in it |
|------|-------------|
| `index.html` | Full 5-page website. Home, Analyzer, Learn, Applications, About |
| `css/style.css` | All colour variables, layout, cards, table styles, animations |
| `js/app.js` | Navigation, demo datasets, analysis pipeline simulation, heatmap drawing (Canvas), SV table, CSV export |
| `python/sv_detector.py` | Full Python pipeline: matrix loading, O/E normalization, feature extraction, Isolation Forest, z-score scan, clustering, classification, report, matplotlib figure |
| `data/sample_matrix_normal.csv` | A small 20×20 matrix you can load to test the Python script |

---

## Technologies used

**Frontend (website)**
- HTML5 + CSS3 + Vanilla JavaScript
- Google Fonts: Syne, JetBrains Mono, Inter
- HTML Canvas API for heatmap drawing
- No frameworks, no build step — just open and run

**Python pipeline**
- `numpy` / `scipy` — matrix math and statistics
- `scikit-learn` — Isolation Forest anomaly detection
- `pandas` — data loading and CSV output
- `matplotlib` — heatmap visualizations

**Real bioinformatics tools (for actual Hi-C processing)**
- Bowtie2 / HiCUP / Juicer — alignment
- cooler / HiC-Pro — contact matrix construction
- Snakemake — workflow management

---

*Project by: B Chandra Deepika | GNITS | B.Tech CSE (Data Science) | 2nd Year*
