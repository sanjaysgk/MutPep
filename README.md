# MutPepGen

Mutation-derived peptide database generator for LC-MS/MS database search. Part of the CAN-IMMUNE project.

## Requirements

- Python 3.10 or higher
- Tkinter (included with most Python installations)

## Installation

```bash
git clone https://github.com/sanjaysgk/MutPepGen.git
cd MutPepGen
python -m venv .venv
source .venv/bin/activate
pip install -e ".[all]"
```

## Usage

```bash
mutpepgen
```

Or run directly:

```bash
python -m mutpepgen.mutpepgen
```

## Input

- Mutation data file (CSV, TSV, or MAF format) with transcript IDs and protein mutations
- Reference protein sequence database (FASTA format)

## Output

- Mutant peptide sequences (FASTA)
- Analysis summary report (HTML)

## License

Monash University
