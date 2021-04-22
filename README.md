# yavalath_mcts

### feature

- Bitboard techniques are used everywhere. It accelerates following functions:
   - enumerate all legal (= not making a line-of-3 pieces) moves
   - enumerate all legal moves that make a check
   - determine whether a player can checkmate in 1 move
   - determine whether a player can checkmate in 3 moves
- MCTS (UCT + random playout)
- during a playout, each player can avoid all the blunders(= lose in 4 moves) and find all the killer-moves(= win in 3 moves).
- This software uses some intel intrinsics, but you can easily change it to use only basic instructions. (#define ONLY_BASIC_INSTRUCTION)
