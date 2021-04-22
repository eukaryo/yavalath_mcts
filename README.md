# yavalath_mcts

### feature

- Bitboard techniques are used everywhere. It accelerates the following functions:
   - enumerate all legal (= not making a line-of-3 pieces) moves
   - enumerate all legal moves that make a check
   - determine whether a player can checkmate in 1 move
   - determine whether a player can checkmate in 3 move
- MCTS (UCT + random playout)
- during playout, each player avoids blunder(= lose in 4 move) and find killer-move (= win in  3 move).
