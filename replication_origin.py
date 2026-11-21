"""
Finding the Origin of Replication (ori) -- interactive toolkit.

Run it with:   python3 replication_origin.py
All of the original functions are kept with the same names and arguments,
so they can still be imported or pasted into course autograders.
"""

# ---------------------------------------------------------------------------
# Sample data: the Vibrio cholerae ori region (540 bp)
# ---------------------------------------------------------------------------
ori = ("ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAG"
       "ATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATAT"
       "CGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGAT"
       "TACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACT"
       "GACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCT"
       "TGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGT"
       "TTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC")

# E. coli: the 500 bp window starting at the minimum skew (position 3,923,620).
# The experimentally verified ori starts 37 nucleotides into this window.
ecoli_window = ("AATGATGATGACGTCAAAAGGATCCGGATAAAACATGGTGATTGCCTCGCATAACGCGGTATGAAAATGGATTGAAGCCC"
                "GGGCCGTGGATTCTACTCAACTTTGTCGGCTTGAGAAAGACCTGGGATCCTGGGTATTAAAAAGAAGATCTATTTATTTAG"
                "AGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCT"
                "GGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCA"
                "AAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGACAGAGTTATCCACAGTAGATCGCACG"
                "ATCTGTATACTTATTTGAGTAAATTAACCCACGATCCCAGCCATTCTTCTGCCGGATCTTCCGGAATGTCGTGATCAAGAA"
                "TGTTGATCTTCAGTG")


# ---------------------------------------------------------------------------
# Input checks
# ---------------------------------------------------------------------------
# DNA is often written in lowercase (many FASTA files are), so every function
# converts its input to uppercase first. Without this, PatternCount("acgacg",
# "ACG") would silently return 0. Arguments that make no sense (an empty
# pattern, k < 1, d < 0) raise a ValueError instead of returning an answer.

def _check_pattern(Pattern):
    if len(Pattern) == 0:
        raise ValueError("Pattern can't be empty.")


def _check_k(k):
    if k < 1:
        raise ValueError("k must be at least 1.")


def _check_d(d):
    if d < 0:
        raise ValueError("d (the number of mismatches) can't be negative.")


# ---------------------------------------------------------------------------
# Part 1: counting k-mers
# ---------------------------------------------------------------------------

# PatternCount() counts how many times Pattern appears in Text (overlaps included).
def PatternCount(Text, Pattern):
    _check_pattern(Pattern)
    Text, Pattern = Text.upper(), Pattern.upper()
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1
    return count


# FrequencyMap() returns a dictionary of every k-mer in Text and how often it appears.
def FrequencyMap(Text, k):
    _check_k(k)
    Text = Text.upper()
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i+k]
        freq[Pattern] = freq.get(Pattern, 0) + 1   # one pass instead of two
    return freq


# FrequentWords() returns the most frequent k-mer(s) in Text.
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    if not freq:            # k is longer than Text -> no k-mers at all
        return words
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


# ---------------------------------------------------------------------------
# Part 2: reverse complements and pattern positions
# ---------------------------------------------------------------------------

def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)      # reverse all letters in a string
    Pattern = Complement(Pattern)   # complement each letter in a string
    return Pattern


def Reverse(Pattern):
    rev = ''
    for char in Pattern:
        rev = char + rev
    return rev


def Complement(Pattern):
    pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complement = ""
    for char in Pattern.upper():
        complement += pairs.get(char, char)   # leave N or other symbols unchanged
    return complement


# PatternMatching() returns every starting position where Pattern appears in Genome.
def PatternMatching(Pattern, Genome):
    _check_pattern(Pattern)
    Pattern, Genome = Pattern.upper(), Genome.upper()
    positions = []
    Pattern_length = len(Pattern)
    Genome_length = len(Genome)
    for i in range(Genome_length - Pattern_length + 1):
        if Genome[i:i + Pattern_length] == Pattern:
            positions.append(i)
    return positions


# ---------------------------------------------------------------------------
# Part 3: where does replication begin?
# ---------------------------------------------------------------------------
# Regular DNA replication shows that:
#   - The forward half-strand is copied as the lagging strand (in Okazaki
#     fragments), so it spends a long time single-stranded and is more prone
#     to mutation.
#   - The reverse half-strand is copied as the leading strand, so it is
#     single-stranded only briefly and is less prone to mutation.
#
# Single-stranded DNA undergoes deamination much faster, where cytosine is
# converted to thymine:
#       C-G -> T-G -> T-A
#
# Over time the forward half-strand loses cytosines, so:
#   FEWER C's in a window == forward half-strand
#   MORE  C's in a window == reverse half-strand

# SymbolArray() counts how many times `symbol` occurs in each window of length
# len(Genome)//2 starting at position i. The genome is treated as circular
# (bacterial chromosomes are circular), which is why it is "extended" by half.
def SymbolArray(Genome, symbol):
    Genome, symbol = Genome.upper(), symbol.upper()
    array = {}
    n = len(Genome)
    if n == 0:
        return array
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(Genome[0:n//2], symbol)
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:              # base leaving the window
            array[i] = array[i] - 1
        if ExtendedGenome[i+(n//2)-1] == symbol:       # base entering the window
            array[i] = array[i] + 1
    return array
# The window with the fewest C's mostly covers the forward half-strand, so the
# start of that window is roughly where the reverse half-strand turns into the
# forward half-strand -- the approximate location of ori.

# Guanine (G) also helps:
#   MORE G's == forward half-strand,  FEWER G's == reverse half-strand
#
# So we track the running total of (#G - #C) as we walk along the genome:
#   skew INCREASES -> forward half-strand
#   skew DECREASES -> reverse half-strand

# SkewArray() returns the running (#G - #C) at every position 0..len(Genome).
def SkewArray(Genome):
    skew = [0]
    for base in Genome.upper():
        if base == "G":
            skew.append(skew[-1] + 1)
        elif base == "C":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

# Plotted as a skew diagram (x = position, y = skew), the minimum point is
# where the skew stops decreasing and starts increasing -- the predicted ori.

# MinimumSkew() returns every position where the skew reaches its minimum.
def MinimumSkew(Genome):
    skew = SkewArray(Genome)
    min_skew = min(skew)
    return [i for i, value in enumerate(skew) if value == min_skew]


# ---------------------------------------------------------------------------
# Part 4: finding DnaA boxes (allowing mismatches)
# ---------------------------------------------------------------------------
# Near the predicted ori you can look for DnaA boxes. The exact 9-mer may not
# appear 3+ times, but counting the 9-mer, its reverse complement, and
# versions with 1 mismatch shows several hits -- DnaA can bind all of these.

# HammingDistance() counts positions where two equal-length strings differ.
def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("Hamming distance needs two strings of the same length.")
    p, q = p.upper(), q.upper()
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance


# ApproximatePatternMatching() returns every position where Pattern appears in
# Text with at most d mismatches.
def ApproximatePatternMatching(Text, Pattern, d):
    _check_pattern(Pattern)
    _check_d(d)
    Text, Pattern = Text.upper(), Pattern.upper()
    positions = []
    pattern_length = len(Pattern)
    text_length = len(Text)
    for i in range(text_length - pattern_length + 1):
        if HammingDistance(Text[i:i+pattern_length], Pattern) <= d:
            positions.append(i)
    return positions


# ApproximatePatternCount() counts how many times Pattern appears in Text with
# at most d mismatches.
def ApproximatePatternCount(Pattern, Text, d):
    _check_pattern(Pattern)
    _check_d(d)
    Text, Pattern = Text.upper(), Pattern.upper()
    count = 0
    pattern_length = len(Pattern)
    text_length = len(Text)
    for i in range(text_length - pattern_length + 1):
        if HammingDistance(Text[i:i+pattern_length], Pattern) <= d:
            count += 1
    return count

# Note: the minimum skew points to a region, not an exact spot -- random
# fluctuations shift it, so always search a window (~500 bp) around it.
# Also, if a genome has a net excess of C over G, the skew ends below 0 and
# the global minimum can fall at the very end of the sequence. Because the
# chromosome is circular, read the diagram as a loop: the ori is the point
# where a long decrease (wrapping around the end) turns into an increase.
# For V. cholerae that turning point is ~153 kb, next to the known ori.


# ---------------------------------------------------------------------------
# Part 4b: most frequent k-mers WITH mismatches (and reverse complements)
# ---------------------------------------------------------------------------
# In E. coli, a 500 bp window starting at the minimum skew has NO 9-mer that
# appears 3+ times exactly (even counting reverse complements). But DnaA also
# binds slightly modified boxes, so we allow d mismatches. The most frequent
# 9-mer with 1 mismatch then turns out to be the real E. coli DnaA box,
# TTATCCACA (and its reverse complement TGTGGATAA).
#
# Important: the most frequent "approximate" k-mer does not have to appear
# exactly in the text at all. So we can't only score k-mers we see -- we score
# every NEIGHBOR (every k-mer within d mismatches) of every k-mer we see.

# Neighbors() returns every string within Hamming distance d of Pattern.
# The course version is recursive (neighbors of Pattern are built from the
# neighbors of Pattern[1:]). That hits Python's recursion limit for patterns
# longer than ~1000 bases and gets very slow for long patterns. This version
# builds each neighbor directly: pick which i positions change (i = 1..d),
# then pick a different base for each of those positions.
def Neighbors(Pattern, d):
    from itertools import combinations, product
    _check_d(d)
    Pattern = Pattern.upper()
    neighborhood = {Pattern}
    for i in range(1, min(d, len(Pattern)) + 1):
        for positions in combinations(range(len(Pattern)), i):
            choices = [[b for b in "ACGT" if b != Pattern[p]] for p in positions]
            for new_bases in product(*choices):
                neighbor = list(Pattern)
                for p, b in zip(positions, new_bases):
                    neighbor[p] = b
                neighborhood.add("".join(neighbor))
    return neighborhood


def NeighborhoodSize(k, d):
    """How many k-mers are within d mismatches of a k-mer (for time estimates)."""
    total, ways = 0, 1                        # ways = C(k, i) * 3^i
    for i in range(min(d, k) + 1):
        total += ways
        ways = ways * (k - i) * 3 // (i + 1)
    return total


# FrequentWordsWithMismatches() returns the k-mer(s) that appear most often in
# Text with at most d mismatches.
def FrequentWordsWithMismatches(Text, k, d):
    _check_k(k)
    _check_d(d)
    Text = Text.upper()
    counts = {}
    for i in range(len(Text) - k + 1):
        for neighbor in Neighbors(Text[i:i+k], d):
            counts[neighbor] = counts.get(neighbor, 0) + 1
    if not counts:
        return []
    m = max(counts.values())
    return sorted(kmer for kmer in counts if counts[kmer] == m)


# Same idea, but each k-mer's score = its count + its reverse complement's
# count, because DnaA boxes occur on both strands.
def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    counts = FrequentWordsWithMismatchesScores(Text, k, d)
    if not counts:
        return []
    m = max(counts.values())
    return sorted(kmer for kmer in counts if counts[kmer] == m)


def FrequentWordsWithMismatchesScores(Text, k, d):
    """Score for every candidate k-mer: Count_d(Pattern) + Count_d(ReverseComplement(Pattern))."""
    _check_k(k)
    _check_d(d)
    Text = Text.upper()
    counts = {}
    for i in range(len(Text) - k + 1):
        window = Text[i:i+k]
        # window is within d of Pattern            -> counts toward Pattern
        for neighbor in Neighbors(window, d):
            counts[neighbor] = counts.get(neighbor, 0) + 1
        # window is within d of ReverseComplement(Pattern), i.e. the reverse
        # complement of window is within d of Pattern -> also counts
        for neighbor in Neighbors(ReverseComplement(window), d):
            counts[neighbor] = counts.get(neighbor, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Part 5: clump finding
# ---------------------------------------------------------------------------
# A k-mer forms an (L, t)-clump if it appears at least t times inside some
# window of length L. DnaA boxes cluster inside the ori, so clumps of 9-mers
# in a ~500 bp window are a strong ori signal -- even when the skew misleads.

# FindClumps() returns every distinct k-mer that forms an (L, t)-clump in Text.
# Instead of rebuilding a FrequencyMap for each window (far too slow for a
# 1.1 million bp genome), it slides the window one base at a time: the k-mer
# leaving the window loses 1 and the k-mer entering it gains 1.
def FindClumps(Text, k, L, t):
    _check_k(k)
    if t < 1:
        raise ValueError("t must be at least 1.")
    Text = Text.upper()
    patterns = set()
    n = len(Text)
    if L > n or k > L:
        return []
    freq = FrequencyMap(Text[0:L], k)        # counts for the first window
    for kmer, count in freq.items():
        if count >= t:
            patterns.add(kmer)
    for i in range(1, n - L + 1):
        leaving = Text[i-1:i-1+k]             # first k-mer of the old window
        freq[leaving] -= 1
        entering = Text[i+L-k:i+L]            # last k-mer of the new window
        freq[entering] = freq.get(entering, 0) + 1
        if freq[entering] >= t:
            patterns.add(entering)
    return sorted(patterns)


# ===========================================================================
# Interactive menu
# ===========================================================================

VALID_BASES = set("ACGTN")


def clean_sequence(raw):
    """Uppercase and strip spaces/newlines/digits; skip FASTA header lines."""
    lines = [line for line in raw.splitlines() if not line.startswith(">")]
    seq = "".join(ch for ch in "".join(lines).upper() if ch.isalpha())
    bad = set(seq) - VALID_BASES
    if bad:
        print(f"  Warning: sequence contains non-DNA letters: {', '.join(sorted(bad))}")
    return seq


def ask_text(prompt):
    while True:
        value = input(prompt).strip()
        if value:
            return value
        print("  Please enter something.")


def ask_dna(prompt):
    while True:
        seq = clean_sequence(input(prompt))
        if seq:
            return seq
        print("  Please enter a DNA sequence (A, C, G, T).")


def ask_int(prompt, minimum=0, maximum=None):
    if maximum is not None and minimum > maximum:
        # No answer could be valid (e.g. k between 1 and 0 for an empty
        # sequence) -- say so instead of asking forever.
        raise ValueError(f"the sequence is too short for this (needs a value from {minimum}, max is {maximum}).")
    while True:
        try:
            value = int(input(prompt))
        except ValueError:
            print("  Please enter a whole number.")
            continue
        if value < minimum or (maximum is not None and value > maximum):
            upper = f" and {maximum}" if maximum is not None else ""
            print(f"  Please enter a number at least {minimum}{upper}.")
            continue
        return value


def show_list(label, items, limit=50):
    print(f"  {label} ({len(items)} total):")
    shown = list(items)[:limit]
    print("   ", " ".join(str(x) for x in shown) if shown else "(none)")
    if len(items) > limit:
        print(f"    ...and {len(items) - limit} more")


# Public copies of the course datasets (from the textbook's website).
VIBRIO_GENOME_URL = "https://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt"
VIBRIO_ORI_URL = "https://bioinformaticsalgorithms.com/data/realdatasets/Replication/v_cholerae_oric.txt"
ECOLI_GENOME_URL = "https://bioinformaticsalgorithms.com/data/realdatasets/Replication/E_coli.txt"


def download_sequence(url):
    """Download a sequence from a link. Returns None if there's no internet."""
    from urllib.request import urlopen
    print(f"  Downloading {url} ...")
    try:
        with urlopen(url, timeout=30) as response:
            seq = clean_sequence(response.read().decode())
    except Exception as e:
        print(f"  Download failed ({e}).")
        return None
    if not seq:
        print("  The download had no sequence in it.")
        return None
    print(f"  Downloaded {len(seq):,} bp.")
    return seq


def load_sequence():
    print("\nChoose a sequence to work with:")
    print("  1) Vibrio cholerae ori (built-in, 540 bp)")
    print("  w) E. coli ori window (built-in, 500 bp from the textbook)")
    print("  2) Vibrio cholerae ori (download from link)")
    print("  3) Vibrio cholerae full genome (download from link, ~1.1 million bp)")
    print("  4) E. coli full genome (download from link, ~4.6 million bp)")
    print("  5) Paste a sequence")
    print("  6) Load from a text or FASTA file")
    choice = input("Choice [1]: ").strip().lower() or "1"
    if choice == "w":
        return ecoli_window
    if choice == "4":
        seq = load_ecoli_genome()
        if seq is None:
            print("  Using the built-in ori instead.")
            return ori
        return seq
    if choice == "2":
        seq = download_sequence(VIBRIO_ORI_URL)
        if seq is None:
            print("  Using the built-in ori instead.")
            return ori
        return seq
    if choice == "3":
        seq = load_vibrio_genome()
        if seq is None:
            print("  Using the built-in ori instead.")
            return ori
        return seq
    if choice == "5":
        return ask_dna("Paste sequence: ")
    if choice == "6":
        seq = read_sequence_file(ask_text("File path: "))
        if seq is None:
            print("  Using the built-in sample instead.")
            return ori
        return seq
    return ori


def read_sequence_file(path):
    """Read a DNA sequence from a text or FASTA file. Returns None (after
    explaining why) if the file can't be opened, isn't text, or is empty."""
    try:
        with open(path) as f:
            seq = clean_sequence(f.read())
    except UnicodeDecodeError:
        print(f"  {path} doesn't look like a text file (is it an image or Word document?).")
        return None
    except OSError as e:
        print(f"  Could not open file: {e}")
        return None
    if not seq:
        print(f"  {path} has no sequence in it.")
        return None
    return seq


def load_vibrio_genome():
    """Load the full V. cholerae genome. In the course this was read in
    'behind the scenes'."""
    return load_dataset("Vibrio_cholerae.txt", VIBRIO_GENOME_URL)


def load_ecoli_genome():
    return load_dataset("E_coli.txt", ECOLI_GENOME_URL)


def load_dataset(filename, url):
    """Use `filename` if it's next to this script; otherwise download it from
    `url` and save a copy so the next run is instant."""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, filename)
    if os.path.exists(path):
        seq = read_sequence_file(path)
        if seq is not None:
            return seq
        print("  Downloading a fresh copy instead.")
    seq = download_sequence(url)
    if seq is None:
        other = input(f"  Enter the path to a local {filename} (or press Enter to skip): ").strip()
        if not other:
            return None
        return read_sequence_file(other)
    try:
        with open(path, "w") as f:
            f.write(seq)
        print(f"  Saved a copy to {path}")
    except OSError:
        pass
    return seq


def estimate_mismatch_search(k, d, windows):
    """Rough (seconds, megabytes) for a mismatch search over `windows` k-mers.
    Measured: ~1 microsecond per neighbor plus ~10 ns per letter, and every
    distinct neighbor is stored (about k + 100 bytes each)."""
    per_kmer = NeighborhoodSize(k, d)
    neighbors = per_kmer * 2 * max(windows, 0)        # x2: reverse complements
    seconds = neighbors * (1e-6 + k * 1e-8)
    distinct = min(neighbors, 4 ** min(k, 30))
    megabytes = distinct * (k + 100) / 1e6
    return seconds, megabytes


def confirm_mismatch_search(k, d, windows):
    """Mismatch searches grow very fast with k and d (k=9, d=3 is ~2,600
    neighbors per k-mer; d=9 is 262,144). Warn before anything slow or
    memory-hungry instead of appearing to hang."""
    seconds, megabytes = estimate_mismatch_search(k, d, windows)
    if seconds < 15 and megabytes < 1000:
        return True
    if seconds < 60:
        time_text = f"about {seconds:.0f} seconds"
    elif seconds < 3600:
        time_text = f"about {seconds / 60:.0f} minute(s)"
    else:
        time_text = f"over {seconds / 3600:,.0f} hour(s)"
    memory_text = f"{megabytes / 1000:,.1f} GB" if megabytes >= 1000 else f"{megabytes:,.0f} MB"
    print(f"  Heads up: k={k}, d={d} means {NeighborhoodSize(k, d):,} neighbors per k-mer.")
    print(f"  This search would take {time_text} and roughly {memory_text} of memory.")
    print("  (DnaA boxes use k=9, d=1, which takes well under a second.)")
    return input("  Run it anyway? [y/N]: ").strip().lower() == "y"


def run_course_examples():
    # E. coli: no DnaA box shows up exactly, but one does with 1 mismatch.
    exact = [kmer for kmer in FrequencyMap(ecoli_window, 9)
             if PatternCount(ecoli_window, kmer)
             + PatternCount(ecoli_window, ReverseComplement(kmer)) >= 3]
    print(f"  E. coli window: 9-mers appearing 3+ times exactly (with reverse complement): {exact or 'none'}")
    best, tied, observed = dnaa_candidates(ecoli_window, 9, 1)
    print(f"  With 1 mismatch, most frequent (score {best}) and present in the window:")
    print(f"  {observed}")
    print("  -> includes TTATCCACA / TGTGGATAA, the experimentally confirmed E. coli DnaA box\n")

    # V. cholerae: approximate DnaA boxes in the ori.
    hits = (ApproximatePatternCount("ATGATCAAG", ori, 1)
            + ApproximatePatternCount("CTTGATCAT", ori, 1))
    print(f"  V. cholerae ori: ATGATCAAG + CTTGATCAT with <= 1 mismatch = {hits} hits")
    print("  (8 binding sites -- the site at position 247 matches both strands' patterns)\n")

    # Exercise 1: count "TGATCA" in the oriC of Vibrio cholerae.
    Text = ori
    Pattern = "TGATCA"
    print(f"  PatternCount(ori, \"{Pattern}\") = {PatternCount(Text, Pattern)}")

    # Exercise 2: find every position of "CTTGATCAT" in the full genome.
    Genome = load_vibrio_genome()
    if Genome is None:
        return
    Pattern = "CTTGATCAT"
    positions = PatternMatching(Pattern, Genome)
    print(f"  PatternMatching(\"{Pattern}\", Genome) =")
    print(f"  {positions}")


def plot_skew(genome):
    skew = SkewArray(genome)
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib isn't installed (pip install matplotlib).")
        print("  Text summary instead:")
        print(f"    min skew = {min(skew)} at {MinimumSkew(genome)}")
        print(f"    max skew = {max(skew)}")
        return
    minima = MinimumSkew(genome)
    plt.figure(figsize=(10, 4))
    # A screen can't show millions of points, so plot at most ~100,000
    # (keeps full-genome plots fast and light on memory).
    step = max(1, len(skew) // 100_000)
    plt.plot(range(0, len(skew), step), skew[::step], linewidth=1)
    plt.scatter(minima, [skew[i] for i in minima], color="red", zorder=3,
                label=f"minimum skew (predicted ori): {minima[:5]}")
    plt.xlabel("Position in genome")
    plt.ylabel("Skew (#G - #C)")
    plt.title("Skew diagram")
    plt.legend()
    plt.tight_layout()
    filename = "skew_diagram.png"
    plt.savefig(filename, dpi=150)
    print(f"  Saved plot to {filename}")
    try:
        plt.show()
    except Exception:
        pass


def dnaa_candidates(region, k, d):
    """Most frequent k-mers (with d mismatches, counting reverse complements)
    in a region. Returns (best score, every tied k-mer, the tied k-mers that
    also appear exactly in the region)."""
    scores = FrequentWordsWithMismatchesScores(region, k, d)
    if not scores:
        return 0, [], []
    best = max(scores.values())
    tied = sorted(kmer for kmer in scores if scores[kmer] == best)
    present = set(FrequencyMap(region, k))
    observed = [kmer for kmer in tied if kmer in present]
    return best, tied, observed


def find_dnaa_boxes(genome):
    """Combined workflow: min skew -> window near it -> frequent 9-mers."""
    window = min(500, len(genome))
    k = ask_int("k-mer length [suggest 9]: ", 1, window)
    d = ask_int("Max mismatches d [suggest 1]: ", 0, k)
    if not confirm_mismatch_search(k, d, window - k + 1):
        print("  Cancelled.")
        return
    min_pos = MinimumSkew(genome)[0]
    print(f"\n  Minimum skew is at position {min_pos:,}.")
    print("  Where should the 500 bp window go?")
    print("    s) start at the minimum skew (what the textbook does for E. coli)")
    print("    c) centred on it")
    print("    e) end at it")
    print("    or type a position to start the window there")
    where = input("  Choice [s]: ").strip().lower() or "s"
    if where == "c":
        start = min_pos - window // 2
    elif where == "e":
        start = min_pos - window
    elif where.lstrip("-").isdigit():
        start = int(where)
    else:
        start = min_pos
    # Bacterial genomes are circular, so let the window wrap past the end.
    start %= len(genome)
    region = (genome + genome)[start:start + window]
    end = (start + window) % len(genome)
    print(f"  Searching positions {start:,}..{end:,} (wraps around the end if needed)")

    best, tied, observed = dnaa_candidates(region, k, d)
    print(f"\n  Most frequent {k}-mers with <= {d} mismatches (+ reverse complements): "
          f"score {best}, {len(tied)} tied")
    print("  Of those, the ones that also appear exactly in the window:")
    shown = set()
    for kmer in observed:
        rc = ReverseComplement(kmer)
        if rc in shown:
            continue
        shown.add(kmer)
        print(f"    {kmer} / {rc}")
    if len(tied) > len(observed):
        show_list("All tied k-mers (including ones never seen exactly)", tied, limit=40)


def clump_report(genome):
    k = ask_int("k-mer length k [suggest 9]: ", 1, len(genome))
    L = ask_int("Window length L [suggest 500]: ", k, len(genome))
    t = ask_int("Minimum times t [suggest 3]: ", 1)
    clumps = FindClumps(genome, k, L, t)
    show_list(f"k-mers forming ({L}, {t})-clumps", clumps)
    # DnaA boxes sit on both strands, so the best ori candidates are
    # k-mers whose reverse complement ALSO forms a clump.
    clump_set = set(clumps)
    paired = sorted({min(c, ReverseComplement(c)) for c in clumps
                     if ReverseComplement(c) in clump_set
                     and len(set(c)) > 1})          # skip runs like CCCCCCCCC
    if not paired:
        print("  No k-mers whose reverse complement also clumps.")
        return
    print("\n  Candidate DnaA boxes (k-mer AND its reverse complement both clump):")
    for kmer in paired:
        rc = ReverseComplement(kmer)
        hits = sorted(PatternMatching(kmer, genome) + PatternMatching(rc, genome))
        print(f"    {kmer} / {rc}: {len(hits)} hits")
        # find the densest window so you can see where the clump is
        best_start, best_count = hits[0], 0
        j = 0
        for i in range(len(hits)):
            while hits[i] - hits[j] > L - k:
                j += 1
            if i - j + 1 > best_count:
                best_count, best_start = i - j + 1, hits[j]
        print(f"      densest cluster: {best_count} hits starting at position {best_start:,}")


def menu():
    genome = load_sequence()
    actions = {
        "1":  "Count a pattern (PatternCount)",
        "2":  "Most frequent k-mers (FrequentWords)",
        "3":  "Reverse complement",
        "4":  "Positions of a pattern (PatternMatching)",
        "5":  "Symbol count per half-genome window (SymbolArray)",
        "6":  "Skew array",
        "7":  "Minimum skew positions (predicted ori)",
        "8":  "Plot skew diagram",
        "9":  "Hamming distance between two strings",
        "10": "Approximate pattern positions (with mismatches)",
        "11": "Approximate pattern count (with mismatches)",
        "12": "Find candidate DnaA boxes near ori (full workflow)",
        "13": "Find clumps (FindClumps) and candidate DnaA boxes",
        "14": "Most frequent k-mers with mismatches (+ reverse complements)",
        "e":  "Run course examples (TGATCA count, CTTGATCAT positions)",
        "s":  "Switch sequence",
        "q":  "Quit",
    }
    while True:
        print(f"\nCurrent sequence: {len(genome):,} bp  ({genome[:30]}{'...' if len(genome) > 30 else ''})")
        for key, label in actions.items():
            print(f"  {key:>2}) {label}")
        choice = input("Choose an option: ").strip().lower()
        print()

        try:
            if choice == "1":
                p = ask_dna("Pattern: ")
                print(f"  {p} appears {PatternCount(genome, p)} times.")
            elif choice == "2":
                k = ask_int("k: ", 1, len(genome))
                words = FrequentWords(genome, k)
                count = PatternCount(genome, words[0]) if words else 0
                show_list(f"Most frequent {k}-mers, each appearing {count} times", words)
            elif choice == "3":
                p = ask_dna("Sequence: ")
                print(f"  Reverse complement: {ReverseComplement(p)}")
            elif choice == "4":
                p = ask_dna("Pattern: ")
                show_list(f"Positions of {p}", PatternMatching(p, genome))
                rc = ReverseComplement(p)
                if rc != p:
                    show_list(f"Positions of reverse complement {rc}", PatternMatching(rc, genome))
            elif choice == "5":
                symbol = ask_text("Symbol (A/C/G/T): ").upper()[0]
                arr = SymbolArray(genome, symbol)
                lowest = min(arr.values())
                lows = [i for i, v in arr.items() if v == lowest]
                print(f"  Min count of {symbol} = {lowest} in window(s) starting at:")
                show_list("Positions", lows)
            elif choice == "6":
                show_list("Skew values", SkewArray(genome))
            elif choice == "7":
                show_list("Minimum skew positions", MinimumSkew(genome))
            elif choice == "8":
                plot_skew(genome)
            elif choice == "9":
                p = ask_dna("String 1: ")
                q = ask_dna("String 2: ")
                print(f"  Hamming distance: {HammingDistance(p, q)}")
            elif choice == "10":
                p = ask_dna("Pattern: ")
                d = ask_int("Max mismatches d: ", 0, len(p))
                show_list(f"Positions of {p} (<= {d} mismatches)",
                          ApproximatePatternMatching(genome, p, d))
            elif choice == "11":
                p = ask_dna("Pattern: ")
                d = ask_int("Max mismatches d: ", 0, len(p))
                print(f"  {p} appears {ApproximatePatternCount(p, genome, d)} times with <= {d} mismatches.")
            elif choice == "12":
                find_dnaa_boxes(genome)
            elif choice == "13":
                clump_report(genome)
            elif choice == "14":
                if len(genome) > 5000:
                    print("  This is meant for a short window (like an ori). Use option 12")
                    print("  to search a 500 bp window of a full genome.")
                    continue
                k = ask_int("k: ", 1, len(genome))
                d = ask_int("Max mismatches d: ", 0, k)
                rc = input("Count reverse complements too? [Y/n]: ").strip().lower() != "n"
                if not confirm_mismatch_search(k, d, len(genome) - k + 1):
                    print("  Cancelled.")
                    continue
                if rc:
                    words = FrequentWordsWithMismatchesAndReverseComplements(genome, k, d)
                else:
                    words = FrequentWordsWithMismatches(genome, k, d)
                show_list(f"Most frequent {k}-mers with <= {d} mismatches", words)
            elif choice == "e":
                run_course_examples()
            elif choice == "s":
                genome = load_sequence()
            elif choice == "q":
                print("Goodbye!")
                break
            else:
                print("  Not a valid option.")
        except ValueError as e:
            print(f"  Error: {e}")
        except KeyboardInterrupt:
            # Ctrl-C during a long search goes back to the menu instead of quitting.
            print("\n  Cancelled.")
        except EOFError:
            raise
        except Exception as e:
            # Anything unexpected: report it and keep the menu running.
            print(f"  Something went wrong: {type(e).__name__}: {e}")


if __name__ == "__main__":
    try:
        menu()
    except (KeyboardInterrupt, EOFError):
        print("\nGoodbye!")
