"""cenfuse: centromeres from satellite arrays, LTR retrotransposons, Hi-C and methylation, cross-validated."""
import os

# numpy asks the kernel for transparent huge pages on large arrays; on a server with fragmented memory every such
# allocation stalls in direct compaction (T2T mammal benchmark: the Hi-C bootstrap ran 30x slower). Results do not
# depend on it. Must be set before numpy is first imported; a value the user set is kept.
os.environ.setdefault("NUMPY_MADVISE_HUGEPAGE", "0")

__version__ = "0.8.0"
