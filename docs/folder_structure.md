# Recommended Folder Structure for Project Documentation

```
docs/
│
├── references/                 # External sources (papers, books, links)
│    ├── papers/                # PDFs of articles (gitignored if large)
│    ├── books/                 # Book chapters, e-books (also gitignored)
│    └── links.md               # Curated list of URLs + comment
│
├── theory/                     # Your learning notes (handwritten knowledge)
│    ├── nbody/                 # Physics, ODEs, integration methods
│    ├── hpc/                   # MPI, OpenMP, load balance, SFCs
│    ├── numerics/              # Floating point, stability, error
│    ├── algorithms/            # Genetic algorithms, PSO, fitness eval
│    └── math/                  # Linear algebra, calculus, etc.
│
├── design/                     # Project architecture
│    ├── architecture.md        # High-level diagram
│    ├── modules/               # Descriptions of each module
│    ├── dataflow/              # Data flow diagrams, pipeline stages
│    ├── decisions/             # ADRs (Architecture Decision Records)
│    └── interface/             # API decisions and rationale
│
├── deliverables/               # Reports, slides, documents to deliver
│    ├── reports/
│    ├── presentations/
│    └── abstract.md
│
│
└── personal-notes/             # Your unfiltered notes, sketches, drafts
     ├── daily-journal.md
     ├── scribbles.md
     └── rough-ideas.md
```