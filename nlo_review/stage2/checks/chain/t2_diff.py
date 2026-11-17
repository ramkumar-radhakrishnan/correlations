# Task 2: Limiting_Cases (L) vs Other_terms (O): character-level diff, row by row
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from tdiff_lc import show
PAIRS = [("GI RowII N1", 41, 41, 45, 45), ("GI RowII N2b", 44, 52, 48, 58), ("GI RowIII N1", 57, 60, 63, 67), ("GI RowIII N2", 73, 78, 70, 76),
         ("GI RowIII N3a", 81, 102, 79, 100), ("GI RowIII N3b", 105, 113, 103, 111), ("GI RowIII A0", 116, 116, 114, 114),
         ("GI RowIV", 137, 138, 118, 121), ("GI RowV 1st", 143, 143, 126, 126), ("GI RowV 2a", 147, 147, 130, 130), ("GI RowV 2b", 150, 159, 133, 142),
         ("GII RowI 1a", 166, 166, 148, 150), ("GII RowI 1b", 169, 169, 153, 155), ("GII RowI 2a", 173, 173, 159, 160), ("GII RowI 2b", 176, 177, 163, 165),
         ("GII RowI 2c", 180, 180, 168, 169), ("GII RowI 2d", 183, 184, 172, 173), ("GII RowI 2e", 187, 188, 176, 177), ("GII RowI 2f", 191, 193, 180, 182),
         ("GII RowIV", 201, 202, 190, 191), ("GIII RowI", 207, 207, 196, 198), ("GIII RowII", 211, 212, 202, 204), ("GIII RowIII 1st", 217, 218, 209, 211),
         ("GIII RowIII 2nd", 222, 227, 215, 222), ("GIII RowV", 233, 234, 228, 230), ("GIII RowVI", 238, 238, 234, 235)]
if __name__ == "__main__":
    for name, a, b, oa, ob in PAIRS:
        print("=====", name); show("L", a, b, "O", oa, ob, ctx=30, maxlen=int(sys.argv[1]) if len(sys.argv) > 1 else 200)
