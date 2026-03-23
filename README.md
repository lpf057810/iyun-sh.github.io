# MR Analysis Pipeline

Unified Mendelian Randomization (MR) analysis pipeline supporting eQTL, pQTL, and Traditional MR analyses.

## Directory Structure

```
MR/
├── scripts/
│   ├── run_mr.sh              # Unified entry script
│   ├── MR_eqtl.R              # eQTL analysis
│   ├── MR_pqtl.R              # pQTL analysis
│   ├── MR_traditional.R       # Traditional MR
│   ├── coloc_from_mr.R        # Colocalization analysis
│   ├── find_common_genes.R    # Find common genes
│   ├── generate_mr_report.R   # Report generator
│   ├── generate_combined_report.R  # Combined eQTL+pQTL report
│   └── config/
│       └── MR.config.example.ini  # Configuration template
├── update_opengwas_token.sh   # Quick token update script
├── logs/                      # Log files (auto-generated)
├── result/                    # Analysis results
│   ├── MR_eqtl/               # eQTL results
│   ├── MR_pqtl/              # pQTL results
│   ├── MR_eqtl_pqtl/         # Combined eQTL+pQTL results
│   ├── MR_traditional/        # Traditional MR results
│   └── coloc/                # Colocalization results
└── report/                   # Word reports (auto-generated)
```

## Usage

### Interactive Mode

```bash
cd /media/desk16/share/secure/MR
./run_mr.sh
```

### Command Line Mode

```bash
# eQTL Analysis
./run_mr.sh -t eqtl -i /path/to/genes.csv -o result/MR_eqtl

# pQTL Analysis with colocalization
./run_mr.sh -t pqtl -i /path/to/genes.csv -o result/MR_pqtl -c 1

# Traditional MR (exposure1 or exposure2 or both)
./run_mr.sh -t traditional -o result/MR_traditional

# eQTL + pQTL Combined Analysis
./run_mr.sh -t eqtl_pqtl -i /path/to/genes.csv -o result/MR_eqtl_pqtl

# eQTL + pQTL + Colocalization
./run_mr.sh -t eqtl_pqtl -i /path/to/genes.csv -o result/MR_coloc -c 1

# With custom parameters
./run_mr.sh -t eqtl -i genes.csv -o result/MR_eqtl -p 1e-6 -k 5000

# With custom token (overrides ~/.Renviron)
./run_mr.sh -t eqtl -i genes.csv -o result/MR_eqtl -T "eyJhbGciOiJSUzI1NiIs..."
```

### Direct R Script Mode

```bash
# eQTL Analysis (direct)
Rscript scripts/MR_eqtl.R -i genes.csv -o results/test_eqtl -p 5e-6 -k 10

# pQTL Analysis (direct)
Rscript scripts/MR_pqtl.R -i genes.csv -o results/test_pqtl -p 5e-6 -k 10

# Traditional MR (direct)
Rscript scripts/MR_traditional.R -o results/test_traditional -p 5e-6 -k 10 -t exposure1
```

## Options

| Option | Description |
|--------|-------------|
| `-t, --type` | MR type: eqtl, pqtl, traditional |
| `-i, --input` | Input gene list file (CSV) |
| `-o, --output` | Output directory |
| `-c, --coloc` | Run colocalization (0 or 1) |
| `-p, --pval` | P-value threshold (default: 5e-8) |
| `-k, --kb` | Distance in kb for clumping (default: 10000) |
| `-T, --token` | OpenGWAS JWT token (optional, overrides env/file) |
| `-h, --help` | Show help message |

## Standardization Features

### 1. Log Files (logs/)
- **Naming**: `MR_{type}.{YYYYMMDD_HHMMSS}.log`
- **Content**: Timestamp, parameters, input/output info, statistics

### 2. Result Files (result/)
- **Subdirectories**: `tables/`, `figures/`
- **Config file**: `*.config.ini` (auto-generated, overwritten each run)
- **Tracking file**: `.{timestamp}` (hidden file for traceability)
- **Plots**: Both `.png` and `.pdf` formats

### 3. Report Files (report/)
- **Naming**: `MR_{type}.{YYYYMMDD_HHMMSS}.docx`
- **Content**: Word report with methods/results/tables/figures

### 4. Configuration
Copy the example config and modify as needed:

```bash
cp scripts/config/MR.config.example.ini scripts/config/MR.config.ini
```

## Token Management

> **Important**: Token update script location: `./update_opengwas_token.sh` (in MR root directory)

### Overview
The pipeline uses OpenGWAS JWT token for API access. Token is loaded with the following priority:

1. **Command line** (highest priority) - Pass via `-T` or `--token` parameter
2. **Environment variable** - `OPENGWAS_JWT` environment variable
3. **Config file** - `~/.Renviron` (lowest priority)

### Quick Token Update

When your token expires, you can update it quickly:

```bash
# Option 1: Direct token string
./update_opengwas_token.sh "eyJhbGciOiJSUzI1NiIs..."

# Option 2: From file
./update_opengwas_token.sh -f /path/to/token.txt

# Option 3: Manual edit
nano ~/.Renviron
# Add: OPENGWAS_JWT=your_token_here
```

### Token Validation

The pipeline automatically validates the token at startup. If invalid, you'll see:

```
[Token Check] FAILED - TOKEN EXPIRED - Please update token
[Token Check] To update token, run: ./update_opengwas_token.sh <new_token>
```

### Retry Mechanism

The pipeline includes automatic retry with exponential backoff for transient API failures:

- **Max retries**: 3 attempts
- **Backoff**: 1s, 2s, 4s (exponential)
- **Token expiry**: Immediate failure (no retry) with clear error message

## Quick Start

### 1. Check Token

```bash
# Run any analysis - token is automatically validated at startup
Rscript scripts/MR_eqtl.R -i test.csv -o results/test
```

If token is invalid, you'll see:
```
[Token Check] FAILED - TOKEN EXPIRED - Please update token
```

### 2. Update Token

```bash
# Quick update
./update_opengwas_token.sh "new_token_here"
```

### 3. Run Analysis

See examples above in "Command Line Mode" section.

## Output Files

Each analysis produces:

| File Type | Description |
|-----------|-------------|
| `00.Complete_MR_Results.csv` | Complete MR results (primary table source) |
| `01.MR_Results_All_Genes.csv` | Per-method MR results for all genes |
| `02.Heterogeneity_All_Genes.csv` | Heterogeneity results |
| `03.Pleiotropy_All_Genes.csv` | Pleiotropy results |
| `06.Three_Filter_Results.csv` | Triple-filter results |
| `07.MR_res_gene.csv` | Final causal gene list |
| `tables/*.csv` | Detailed results per gene/protein |
| `figures/*.png` | Visualization plots |
| `figures/*.pdf` | Publication-ready plots |
| `*_node.qs2` | Node snapshots for reproducibility |
| `logs/*sessionInfo*.txt` | Session info snapshots |

## Requirements

- R (>= 4.0)
- Required R packages:
  - TwoSampleMR
  - ieugwasr
  - ggplot2
  - dplyr
  - data.table
  - MRPRESSO
  - ggsci
  - optparse
