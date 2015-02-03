# sample generic secondary structure prediction file read by seq2trj
#
# comments must begin with '#'
#
# columns should be space separated, and are in the following order:
#
# 1-letter amino acid code
# 1-state secondary structure prediction ("H", "E" or "C")
# pH (probability of helix)
# pE (probability of extended/strand)
# pC (probability of coil)
#
# pH+pE+pC=100 should hold for all rows and pH, pE and pC should all
# be integers (or floats) in [0,100]
#
# all residues in the protein should be listed, in order, one residue
# per row
#
# if your secondary structure prediction method provides only a one-
# state prediction (no columns 3, 4 and 5) then when creating this
# generic secondary structure file just fill in 0, 0, 100 for the
# probabilities (with the 100 corrsponding to the 1-state prediction
# of course)
#
M C 1 1 98
V C 12 20 68
L C 13 15 72
S C 16 17 67
E H 49 15 36
G H 58 12 30
E H 56 16 28
W H 56 29 15
Q H 68 28 4
L H 76 23 1
V H 68 29 3
L H 67 28 5
H H 76 19 5
V H 73 18 9
W H 71 17 12
A H 69 21 10
K H 76 14 10
V H 81 15 4
E H 81 15 4
A H 82 12 6
D H 87 2 11
V H 92 1 7
A H 87 0 13
G H 83 0 17
H H 73 0 27
G H 91 0 9
Q H 92 3 5
D H 93 5 2
I H 97 2 1
L H 97 3 0
I H 97 2 1
R H 96 3 1
L H 96 2 2
F H 92 2 6
K H 84 1 15
S H 71 0 29
H C 28 0 72
P H 70 0 30
E H 90 1 9
T H 95 1 4
L H 96 1 3
E H 97 1 2
K H 91 1 8
F H 74 1 25
D C 48 0 52
R H 82 1 17
F H 84 3 13
K H 87 3 10
H H 83 1 16
L H 83 0 17
K H 91 1 8
T H 93 1 6
E H 95 1 4
A H 96 1 3
E H 97 1 2
M H 95 0 5
K H 96 1 3
A H 97 1 2
S H 98 1 1
E H 98 1 1
D H 96 0 4
L H 88 0 12
K H 74 0 26
K H 66 0 34
H C 29 0 71
G C 37 0 63
V H 55 3 42
T H 67 21 12
V H 80 14 6
L H 95 5 0
T H 97 1 2
A H 96 1 3
L H 94 2 4
G H 95 2 3
A H 97 2 1
I H 98 1 1
L H 98 1 1
K H 96 0 4
K H 83 1 16
K H 51 0 49
G C 16 0 84
H C 10 1 89
H C 27 2 71
E H 65 2 33
A H 79 4 17
E H 84 4 12
L H 86 1 13
K H 83 1 16
P H 94 1 5
L H 96 1 3
A H 97 1 2
Q H 97 1 2
S H 94 1 5
H H 91 1 8
A H 68 0 32
T H 67 0 33
K C 43 0 57
H C 24 1 75
K C 15 4 81
I C 12 4 84
P C 40 7 53
I H 85 4 11
K H 85 3 12
Y H 94 1 5
L H 96 1 3
E H 97 0 3
F H 98 1 1
I H 95 1 4
S H 96 1 3
E H 96 1 3
A H 96 1 3
I H 97 1 2
I H 96 2 2
H H 87 9 4
V H 90 6 4
L H 91 6 3
H H 76 16 8
S H 69 13 18
R H 59 6 35
H C 20 2 78
P C 35 1 64
G C 36 1 63
N C 48 3 49
F C 20 3 77
G C 27 4 69
A H 79 1 20
D H 83 1 16
A H 93 2 5
Q H 91 1 8
G H 96 1 3
A H 96 1 3
M H 97 1 2
N H 97 1 2
K H 97 1 2
A H 97 1 2
L H 98 0 2
E H 97 1 2
L H 96 1 3
F H 94 1 5
R H 95 1 4
K H 95 1 4
D H 96 1 3
I H 97 1 2
A H 98 1 1
A H 98 1 1
K H 98 1 1
Y H 98 1 1
K H 96 1 3
E H 94 2 4
L H 81 1 18
G H 53 1 46
Y C 37 3 60
Q C 14 4 82
G C 2 0 98
