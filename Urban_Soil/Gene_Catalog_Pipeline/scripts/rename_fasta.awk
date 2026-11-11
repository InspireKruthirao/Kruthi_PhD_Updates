# Rename FASTA headers to <p>.000_000_000, ... in file order.
# -v p=PREFIX   new ID prefix
# -v keep=1     keep the old ID after the new one (">new old")
# -v ids=FILE   old IDs are written here, one per line, in order
/^>/ {
    old = substr($1, 2)
    print old > ids
    printf ">%s.%03d_%03d_%03d%s\n", p, int(n / 1000000), int(n / 1000) % 1000, n % 1000, (keep ? " " old : "")
    n++
    next
}
{ print }
