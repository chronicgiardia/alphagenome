# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

AlphaGenome is a Python SDK for interacting with Google DeepMind's AlphaGenome API, a DNA sequence model that provides multimodal predictions for genomic analysis. The SDK enables:
- DNA sequence and genomic interval predictions
- Variant effect prediction and scoring
- Visualization of genomic predictions (RNA-seq, ATAC-seq, ChIP-seq, contact maps, etc.)
- Support for human (Homo sapiens) and mouse (Mus musculus) genomes

## Build and Development Commands

### Installation
```bash
# Standard installation from local clone
pip install ./alphagenome

# Development installation with optional dependencies
pip install ./alphagenome[dev]
pip install ./alphagenome[docs]
```

### Testing
```bash
# Run all tests using hatch
python -m hatch test

# Run tests for specific Python version
python -m hatch test --python 3.10
python -m hatch test --python 3.13

# Run tests on all supported Python versions (3.10-3.13)
python -m hatch test --all
```

### Code Quality Checks
```bash
# Run all checks (formatting and linting)
python -m hatch run check:all

# Run formatting check only (pyink - Google style formatter)
python -m hatch run check:format

# Run linting only (pylint)
python -m hatch run check:lint
```

### Documentation
```bash
# Build documentation using Sphinx
cd docs
make html

# Documentation will be in docs/build/html/
```

### Build
The project uses `hatch` with a custom build hook (`hatch_build.py`) that auto-generates Python bindings from protobuf definitions during the build process. The proto files are in `src/alphagenome/protos/`.

## Architecture and Code Structure

### Core Package Structure
- `src/alphagenome/` - Main package
  - `data/` - Data structures and utilities
    - `genome.py` - Genomic intervals, variants, and coordinate systems
    - `track_data.py` - 1D genomic track data (e.g., RNA-seq, ATAC-seq)
    - `junction_data.py` - Splice junction data structures
    - `gene_annotation.py` - Gene and transcript annotations
    - `ontology.py` - Biological ontology term handling
    - `fold_intervals.py` - Genome fold definitions for cross-validation
  - `models/` - Model interaction and prediction
    - `dna_client.py` - Main client for AlphaGenome API (gRPC)
    - `dna_model.py` - Abstract base class defining the model interface
    - `dna_output.py` - Output types and data structures
    - `variant_scorers.py` - Variant effect scoring strategies
    - `interval_scorers.py` - Genomic interval scoring
  - `visualization/` - Plotting and visualization
    - `plot_components.py` - Composable plotting components
    - `plot_transcripts.py` - Transcript and gene visualization
    - `plot.py` - High-level plotting utilities
  - `interpretation/` - Model interpretation tools
    - `ism.py` - In silico mutagenesis (ISM) analysis
  - `protos/` - Auto-generated gRPC/protobuf code (DO NOT EDIT)

### Key Design Patterns

**API Client Architecture:**
- `dna_client.py` implements gRPC streaming client with retry logic
- Uses `retry_rpc` decorator for automatic retries on `RESOURCE_EXHAUSTED` and `UNAVAILABLE` errors
- Supports multiple model versions (ALL_FOLDS, FOLD_0-3) for cross-validation
- Maximum 20 variant scorers per request (`MAX_VARIANT_SCORERS_PER_REQUEST`)

**Genomic Intervals:**
- All intervals use 0-based half-open coordinates `[start, end)`
- Strand notation: `+` (forward), `-` (reverse), `.` (unstranded)
- Intervals support operations: `resize()`, `center()`, `clip()`, `split()`, `expand()`
- String format: `chr22:35677410-36725986:+`

**Variant Representation:**
- Variants use 0-based positions
- Reference and alternate bases must be from `{A, C, G, T, N}`
- Supports SNPs, insertions, deletions, and complex variants

**Prediction Outputs:**
- Multiple output modalities: `RNA_SEQ`, `ATAC`, `DNASE`, `CAGE`, `PROCAP`, `CHIP_HISTONE`, `CHIP_TF`, `CONTACT_MAPS`, `SPLICE_SITES`, `SPLICE_JUNCTION_COUNTS`
- Track data stored as numpy arrays with associated interval metadata
- Ontology terms (e.g., `UBERON:0001157`) specify tissue/cell type context

**Variant Scoring Strategies:**
- `CenterMaskScorer` - Aggregates predictions in a window around the variant
- `ContactMapScorer` - Handles 2D contact map predictions
- `GeneMaskScorer` - Gene-level aggregation with different modes (LFC, ACTIVE, SPLICING)
- `SpliceJunctionScorer` - Evaluates splice junction effects
- Aggregation types: `DIFF_MEAN`, `DIFF_SUM`, `L2_DIFF`, etc.

### Testing Conventions
- Test files follow `*_test.py` naming convention
- Tests use `absltest` framework
- `conftest.py` configures FLAGS for test environment
- Matplotlib backend set to `'agg'` in test environment via `MPLBACKEND` env var

## Code Style and Standards

### Google Python Style Guide
- 80 character line length
- 2-space indentation (enforced by pyink)
- Use majority quotes (enforced by pyink)
- Function/variable names: `snake_case`
- Class names: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Pylint configuration in `.pylintrc`

### Type Annotations
- Python 3.10+ required
- Use type hints for function signatures
- `jaxtyping` used for array shape annotations
- Import from `typing_extensions` for backward compatibility features

### Protobuf Code
- Never manually edit files in `src/alphagenome/protos/`
- Auto-generated during build from `.proto` files
- Regenerate by triggering build process

## Common Workflows

### Making Predictions
```python
from alphagenome.data import genome
from alphagenome.models import dna_client

model = dna_client.create(API_KEY)

# Predict on sequence
output = model.predict_sequence(
    sequence="ACGT...",
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157']
)

# Predict on genomic interval
interval = genome.Interval(chromosome='chr22', start=100000, end=116384)
output = model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157']
)

# Predict variant effect
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C'
)
variant_output = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157']
)
```

### Variant Scoring
```python
from alphagenome.models import variant_scorers

# Define scorer
scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.RNA_SEQ,
    width=10001,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN
)

# Score variants
scores = model.score_variants(
    interval=interval,
    variants=[variant1, variant2],
    variant_scorers=[scorer],
    ontology_terms=['UBERON:0001157']
)
```

### Visualization
```python
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt

plot_components.plot(
    [
        plot_components.OverlaidTracks(
            tdata={'REF': ref_track, 'ALT': alt_track},
            colors={'REF': 'dimgrey', 'ALT': 'red'}
        )
    ],
    interval=output_interval,
    annotations=[plot_components.VariantAnnotation([variant])]
)
plt.show()
```

## Important Constraints

### Sequence Length Support
Supported sequence lengths are powers of 2:
- 2KB (2,048 bp) - `SEQUENCE_LENGTH_2KB`
- 16KB (16,384 bp) - `SEQUENCE_LENGTH_16KB`
- 100KB (131,072 bp) - `SEQUENCE_LENGTH_100KB`
- 500KB (524,288 bp) - `SEQUENCE_LENGTH_500KB`
- 1MB (1,048,576 bp) - `SEQUENCE_LENGTH_1MB`

### API Limitations
- ISM intervals automatically chunked to max width of 10 bp
- Maximum 20 variant scorers per request
- API designed for 1000s-100000s of predictions, not millions
- Rate limits vary based on demand

### Organism Support
- Most features support human (`HOMO_SAPIENS`) and mouse (`MUS_MUSCULUS`)
- PA_QTL scoring only available for human

## References
- Documentation: https://www.alphagenomedocs.com/
- Community: https://www.alphagenomecommunity.com
- Colab notebooks in `colabs/` directory for examples
- Paper: Avsec et al. 2025 (https://doi.org/10.1101/2025.06.25.661532)
