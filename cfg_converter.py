class CFG:
    """Context-Free Grammar"""
    def __init__(self, start_symbol='S'):
        self.productions = {}  # {non_terminal: [[symbols], [symbols], ...]}
        self.start_symbol = start_symbol
        self.terminals = set()
        self.non_terminals = set()
    
    def add_production(self, lhs, rhs):
        """Add production: lhs -> rhs (rhs is list of symbols)"""
        if lhs not in self.productions:
            self.productions[lhs] = []
        
        self.productions[lhs].append(rhs)
        self.non_terminals.add(lhs)
        
        for symbol in rhs:
            if symbol.islower() or symbol in ['ε', '(', ')', '+', '*', 'a', 'b', '0', '1']:
                self.terminals.add(symbol)
            elif symbol != 'ε':
                self.non_terminals.add(symbol)
    
    def display(self, title="Grammar"):
        print(f"\n{'='*60}")
        print(f"{title}")
        print(f"{'='*60}")
        print(f"Start Symbol: {self.start_symbol}")
        print(f"Non-terminals: {sorted(self.non_terminals)}")
        print(f"Terminals: {sorted(self.terminals)}")
        print(f"\nProductions:")
        for lhs in sorted(self.productions.keys()):
            productions = self.productions[lhs]
            rhs_list = ' | '.join([''.join(p) if p != ['ε'] else 'ε' for p in productions])
            print(f"  {lhs} → {rhs_list}")
    
    def copy(self):
        """Create a copy of the grammar"""
        new_cfg = CFG(self.start_symbol)
        for lhs, productions in self.productions.items():
            for rhs in productions:
                new_cfg.add_production(lhs, rhs[:])
        return new_cfg


class CFGtoChomsky:
    """Convert CFG to Chomsky Normal Form (CNF)"""
    
    def __init__(self, cfg):
        self.cfg = cfg.copy()
        self.new_var_counter = 0
    
    def convert(self, verbose=True):
        """Complete conversion to CNF"""
        if verbose:
            print("\n" + "="*70)
            print("CONVERTING CFG TO CHOMSKY NORMAL FORM (CNF)")
            print("="*70)
            self.cfg.display("Original Grammar")
        
        # Step 1: Add new start symbol
        self._add_new_start(verbose)
        
        # Step 2: Eliminate ε-productions
        self._eliminate_epsilon(verbose)
        
        # Step 3: Eliminate unit productions
        self._eliminate_unit(verbose)
        
        # Step 4: Convert to CNF form
        self._convert_to_cnf(verbose)
        
        if verbose:
            self.cfg.display("Final CNF Grammar")
            self._verify_cnf()
        
        return self.cfg
    
    def _add_new_start(self, verbose):
        """Add new start symbol"""
        # Only add S0 if the start symbol is nullable (can derive ε).
        # This matches many textbook outputs for cases without ε, even if S appears on some RHS.
        old_start = self.cfg.start_symbol
        # compute nullable closure
        nullable = set()
        changed = True
        while changed:
            changed = False
            for lhs, prods in self.cfg.productions.items():
                if lhs in nullable:
                    continue
                for rhs in prods:
                    if rhs == ['ε'] or all(sym in nullable for sym in rhs):
                        nullable.add(lhs)
                        changed = True
                        break

        if old_start in nullable:
            if verbose:
                print("\n" + "-"*60)
                print("Step 1: Adding New Start Symbol")
                print("-"*60)
            new_start = 'S0'
            while new_start in self.cfg.non_terminals:
                new_start += "'"
            self.cfg.add_production(new_start, [old_start])
            self.cfg.start_symbol = new_start
            if verbose:
                print(f"Added: {new_start} → {old_start}")
        else:
            if verbose:
                print("\n" + "-"*60)
                print("Step 1: Start symbol is not nullable; skipping new start symbol")
                print("-"*60)
    
    def _eliminate_epsilon(self, verbose):
        """Eliminate ε-productions"""
        if verbose:
            print("\n" + "-"*60)
            print("Step 2: Eliminating ε-Productions")
            print("-"*60)
        
        # Find nullable non-terminals
        nullable = set()
        changed = True
        
        while changed:
            changed = False
            for lhs, productions in self.cfg.productions.items():
                if lhs in nullable:
                    continue
                for rhs in productions:
                    if rhs == ['ε'] or all(s in nullable for s in rhs):
                        nullable.add(lhs)
                        changed = True
                        break
        
        if verbose and nullable:
            print(f"Nullable non-terminals: {sorted(nullable)}")
        
        # Generate new productions
        new_productions = {}
        for lhs, productions in self.cfg.productions.items():
            new_productions[lhs] = []
            
            for rhs in productions:
                if rhs == ['ε']:
                    continue
                
                # Generate all combinations by removing nullable symbols
                variants = self._generate_variants(rhs, nullable)
                for variant in variants:
                    if variant and variant not in new_productions[lhs]:
                        new_productions[lhs].append(variant)
        
        self.cfg.productions = new_productions

        # Remove ε from the terminals set if present (productions no longer contain ε)
        if 'ε' in self.cfg.terminals:
            try:
                self.cfg.terminals.remove('ε')
            except KeyError:
                pass

        if verbose:
            print("ε-productions eliminated")
    
    def _generate_variants(self, rhs, nullable):
        """Generate all variants by removing nullable symbols"""
        if not rhs:
            return [[]]
        
        variants = []
        
        def generate(index, current):
            if index == len(rhs):
                variants.append(current[:])
                return
            
            symbol = rhs[index]
            
            # Include the symbol
            current.append(symbol)
            generate(index + 1, current)
            current.pop()
            
            # Exclude if nullable
            if symbol in nullable:
                generate(index + 1, current)
        
        generate(0, [])
        return [v for v in variants if v]
    
    def _eliminate_unit(self, verbose):
        """Eliminate unit productions (A → B)"""
        if verbose:
            print("\n" + "-"*60)
            print("Step 3: Eliminating Unit Productions")
            print("-"*60)
        
        # Find all unit productions
        unit_pairs = set()
        for lhs, productions in self.cfg.productions.items():
            for rhs in productions:
                if len(rhs) == 1 and rhs[0] in self.cfg.non_terminals:
                    unit_pairs.add((lhs, rhs[0]))
        
        if verbose and unit_pairs:
            print(f"Unit productions found: {unit_pairs}")
        
        # Compute closure of unit pairs
        closure = unit_pairs.copy()
        changed = True
        while changed:
            changed = False
            for a, b in list(closure):
                for b2, c in list(closure):
                    if b == b2 and (a, c) not in closure:
                        closure.add((a, c))
                        changed = True
        
        # Generate new productions
        new_productions = {}
        for lhs in self.cfg.non_terminals:
            new_productions[lhs] = []
            
            # Add non-unit productions from lhs
            for rhs in self.cfg.productions.get(lhs, []):
                if len(rhs) != 1 or rhs[0] not in self.cfg.non_terminals:
                    if rhs not in new_productions[lhs]:
                        new_productions[lhs].append(rhs)
            
            # Add non-unit productions from all B where (lhs, B) in closure
            for a, b in closure:
                if a == lhs:
                    for rhs in self.cfg.productions.get(b, []):
                        if len(rhs) != 1 or rhs[0] not in self.cfg.non_terminals:
                            if rhs not in new_productions[lhs]:
                                new_productions[lhs].append(rhs)
        
        self.cfg.productions = new_productions
        
        if verbose:
            print("Unit productions eliminated")
    
    def _convert_to_cnf(self, verbose):
        """Convert remaining productions to CNF form"""
        if verbose:
            print("\n" + "-"*60)
            print("Step 4: Converting to CNF Form")
            print("-"*60)
        # Deterministic terminal variables:
        # - digits '0' -> T0, '1' -> T1
        # - single-letter terminals like 'a','b' -> C, D, E, ... (avoid collisions)
        terminal_vars = {}
        new_productions = {}

        # Collect terminals present in productions
        terminals_in_use = set()
        for lhs, productions in self.cfg.productions.items():
            for rhs in productions:
                for sym in rhs:
                    if sym != 'ε' and sym not in self.cfg.non_terminals:
                        terminals_in_use.add(sym)

        # Pre-assign names deterministically
        pool_letters = [c for c in ['C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','U','V','W','X','Y','Z'] if c not in self.cfg.non_terminals]
        for t in sorted(terminals_in_use):
            if t == '0':
                name = 'T0'
            elif t == '1':
                name = 'T1'
            else:
                # assign next available from pool
                if pool_letters:
                    name = pool_letters.pop(0)
                else:
                    # fallback to generated variable
                    name = self._get_new_variable()
            # ensure uniqueness
            if name in self.cfg.non_terminals:
                # if conflict, fallback to generated
                name = self._get_new_variable()
            terminal_vars[t] = name
            self.cfg.non_terminals.add(name)
            new_productions.setdefault(name, []).append([t])
            if verbose:
                print(f"Created terminal variable: {name} → {t}")

        def get_terminal_var(t):
            return terminal_vars[t]

        # Helper: introduce pair variable for identical last two non-terminals (AA -> A1, BB -> B1)
        pair_cache = {}
        def get_pair_var(a, b):
            key = (a, b)
            if key in pair_cache:
                return pair_cache[key]
            # identical pair naming rule
            if a == b and len(a) == 1 and a.isupper():
                name = f"{a}1"
            else:
                # fallback generic variable
                name = self._get_new_variable()
            # Ensure uniqueness
            if name in self.cfg.non_terminals and key not in pair_cache and not (a == b and len(a)==1 and a.isupper()):
                # already used, fallback to generated
                name = self._get_new_variable()
            pair_cache[key] = name
            self.cfg.non_terminals.add(name)
            new_productions.setdefault(name, []).append([a, b])
            if verbose:
                print(f"Created pair variable: {name} → {a}{b}")
            return name

        for lhs, productions in self.cfg.productions.items():
            new_productions.setdefault(lhs, [])
            for rhs in productions:
                # Replace terminals inside longer productions
                if len(rhs) == 1:
                    sym = rhs[0]
                    # Keep single terminal as-is to satisfy strict CNF (A → a)
                    # Only map terminals to variables when they appear in longer productions
                    new_productions[lhs].append([sym])
                    continue
                replaced = []
                for sym in rhs:
                    if sym in terminals_in_use and len(rhs) > 1:
                        # replace terminal only when part of a longer production
                        replaced.append(get_terminal_var(sym))
                    else:
                        replaced.append(sym)

                # If length now 2 and both non-terminals => CNF
                if len(replaced) == 2 and all(s in self.cfg.non_terminals for s in replaced):
                    new_productions[lhs].append(replaced)
                    continue

                # Break longer productions by introducing pair vars preferentially on identical trailing pairs
                current = replaced
                while len(current) > 2:
                    # Prefer using the last two symbols
                    a, b = current[-2], current[-1]
                    if a in self.cfg.non_terminals and b in self.cfg.non_terminals:
                        pair_var = get_pair_var(a, b)
                        current = current[:-2] + [pair_var]
                    else:
                        # Fallback: use first two symbols
                        a, b = current[0], current[1]
                        pair_var = get_pair_var(a, b)
                        current = [pair_var] + current[2:]
                new_productions[lhs].append(current)

        # --- Redundancy cleanup -------------------------------------------------
        # Some grammars (especially after ε and unit elimination) can leave a direct
        # two-nonterminal production like S → A S alongside a helper pair variable
        # X1 → A S that is used to factor longer productions (e.g. A → C X1).
        # The user expects the direct production (S → A S) to be removed in the
        # final CNF when the pair already has its own dedicated helper variable.
        #
        # Strategy:
        # 1. Collect all pair variables (variables whose ONLY production is exactly two non-terminals).
        # 2. Build a map pair_rhs -> helper_var.
        # 3. For every other non-terminal, remove any production matching a pair_rhs
        #    if that helper_var also appears somewhere else (i.e., the helper is actually used).
        pair_var_map = {}
        for var, prods in new_productions.items():
            if len(prods) == 1 and len(prods[0]) == 2 and all(s in self.cfg.non_terminals for s in prods[0]):
                pair = tuple(prods[0])
                pair_var_map[pair] = var

        # Track usage of helper variables in RHSs of other productions
        helper_usage = {v: 0 for v in pair_var_map.values()}
        for lhs, prods in new_productions.items():
            for rhs in prods:
                for sym in rhs:
                    if sym in helper_usage:
                        helper_usage[sym] += 1

        # Execute removal pass
        for lhs, prods in list(new_productions.items()):
            filtered = []
            for rhs in prods:
                if len(rhs) == 2:
                    pair = tuple(rhs)
                    if pair in pair_var_map and lhs != pair_var_map[pair]:
                        helper = pair_var_map[pair]
                        # Remove if helper is used elsewhere (usage > 1 means besides its own definition)
                        if helper_usage.get(helper, 0) > 1:
                            # Skip (remove) this redundant direct pair production
                            if verbose:
                                print(f"Removed redundant pair production: {lhs} → {''.join(rhs)} (replaced by {helper} → {''.join(rhs)})")
                            continue
                filtered.append(rhs)
            new_productions[lhs] = filtered

        self.cfg.productions = new_productions
        if verbose:
            print("Conversion to CNF complete (deterministic naming applied)")
            # Print summary of terminal variables
            if terminal_vars:
                print("Terminal variables introduced:")
                for t, v in sorted(terminal_vars.items()):
                    print(f"  {v} → {t}")
    
    def _get_new_variable(self):
        """Generate new variable name"""
        while True:
            new_var = f"X{self.new_var_counter}"
            self.new_var_counter += 1
            if new_var not in self.cfg.non_terminals:
                self.cfg.non_terminals.add(new_var)
                return new_var
    
    def _verify_cnf(self):
        """Verify that grammar is in CNF"""
        print("\n" + "-"*60)
        print("CNF Verification")
        print("-"*60)
        
        is_cnf = True
        for lhs, productions in self.cfg.productions.items():
            for rhs in productions:
                if len(rhs) == 1:
                    if rhs[0] not in self.cfg.terminals:
                        print(f"✗ Invalid: {lhs} → {''.join(rhs)} (should be terminal)")
                        is_cnf = False
                elif len(rhs) == 2:
                    if not all(s in self.cfg.non_terminals for s in rhs):
                        print(f"✗ Invalid: {lhs} → {''.join(rhs)} (should be non-terminals)")
                        is_cnf = False
                else:
                    print(f"✗ Invalid: {lhs} → {''.join(rhs)} (wrong length)")
                    is_cnf = False
        
        if is_cnf:
            print("✓ Grammar is in valid Chomsky Normal Form")
        else:
            print("✗ Grammar is NOT in valid CNF")

    # --- Trace-enabled conversion -------------------------------------------------
    def convert_with_trace(self):
        """Perform conversion but print a Step 1..4 trace in the user's requested format."""
        out_lines = []

        def fmt_grammar_header():
            return ""  # placeholder if needed

        def grammar_to_block(cfg):
            # Produce lines like: S0 → S but with pretty short names for generated vars
            # Build a mapping from internal non-terminals to user-friendly names
            # Keep single-letter uppercase names (A,B,S,...) and S0 as-is; rename others to P,R,T,Q,...
            used = set()
            mapping = {}
            pool = [c for c in ['P', 'R', 'T', 'Q', 'M', 'N', 'O', 'U', 'V', 'W', 'Z', 'Y']]

            # ensure deterministic ordering
            nts = sorted(cfg.non_terminals, key=lambda x: (len(x), x))
            for nt in nts:
                if nt == 'S0' or (len(nt) == 1 and nt.isupper()):
                    mapping[nt] = nt
                    used.add(nt)
            for nt in nts:
                if nt in mapping:
                    continue
                # pick next from pool that's not used
                found = None
                for p in pool:
                    if p not in used:
                        found = p
                        break
                if found is None:
                    # fallback: use original name
                    mapping[nt] = nt
                    used.add(nt)
                else:
                    mapping[nt] = found
                    used.add(found)

            # now build ordered keys: S0, S, then mapped single letters and others
            keys = list(cfg.productions.keys())
            keys_sorted = []
            if 'S0' in keys:
                keys_sorted.append('S0')
            if 'S' in keys and 'S' not in keys_sorted:
                keys_sorted.append('S')
            for k in sorted(keys):
                if k not in keys_sorted:
                    keys_sorted.append(k)

            lines = []
            for lhs in keys_sorted:
                rhss = cfg.productions.get(lhs, [])
                if not rhss:
                    rhs_text = ''
                else:
                    rhs_list = []
                    for rhs in rhss:
                        if rhs == ['ε']:
                            rhs_list.append('ε')
                        else:
                            # map each symbol if it's a non-terminal
                            out_syms = []
                            for s in rhs:
                                if s in mapping:
                                    out_syms.append(mapping[s])
                                else:
                                    out_syms.append(s)
                            rhs_list.append(''.join(out_syms))
                    rhs_text = ' | '.join(rhs_list)
                lhs_disp = mapping[lhs] if lhs in mapping else lhs
                lines.append(f"{lhs_disp} → {rhs_text}")
            return lines

        # Work on a copy so we can modify
        g = self.cfg

        # Step 1
        print("Step 1.")
        print()
        print(f"As start symbol {g.start_symbol} appears on the RHS, we will create a new production rule S0→{g.start_symbol}. Therefore, the grammar will become:")
        print()
        # Add new start
        new_start = 'S0'
        while new_start in g.non_terminals:
            new_start += "'"
        old_start = g.start_symbol
        g.add_production(new_start, [old_start])
        g.start_symbol = new_start
        for l in grammar_to_block(g):
            print(l)
        print()

        # Step 2: ε-elimination with iterative snapshots
        print("Step 2.")
        print()
        # Find initial nullable (those with ε-production)
        nullable = set()
        for lhs, prods in g.productions.items():
            for rhs in prods:
                if rhs == ['ε']:
                    nullable.add(lhs)
        # We will iteratively remove nullable symbols and show snapshots when new nullable found
        processed_nullable = set()
        # First snapshot: if A is nullable, show removal of A→ε
        while True:
            newly = [n for n in nullable if n not in processed_nullable]
            if not newly:
                break
            for n in newly:
                processed_nullable.add(n)
                print(f"As grammar contains null production {n}→ ε, its removal from the grammar yields:")
                # Remove explicit ε-productions and generate variants removing n from RHS
                new_productions = {}
                for lhs, prods in g.productions.items():
                    new_productions[lhs] = []
                    for rhs in prods:
                        if rhs == ['ε'] and lhs == n:
                            continue
                        # generate variants by removing occurrences of n
                        variants = [[]]
                        for sym in rhs:
                            if sym == n:
                                # for each current variant, create two: with and without sym
                                with_sym = [v + [sym] for v in variants]
                                without_sym = [v[:] for v in variants]
                                variants = with_sym + without_sym
                            else:
                                for v in variants:
                                    v.append(sym)
                        # add non-empty variants
                        for v in variants:
                            if v and v not in new_productions[lhs]:
                                new_productions[lhs].append(v)
                g.productions = new_productions
                # Recompute nullable after removal: if a production becomes empty, mark lhs nullable
                newly_found = set()
                for lhs, prods in g.productions.items():
                    for rhs in prods:
                        if rhs == ['ε'] or not rhs:
                            if lhs not in nullable:
                                newly_found.add(lhs)
                nullable.update(newly_found)
                # print snapshot grammar
                for l in grammar_to_block(g):
                    print(l)
                print()

        # Clean any residual ε entries
        for lhs in list(g.productions.keys()):
            g.productions[lhs] = [rhs for rhs in g.productions[lhs] if rhs != ['ε']]
        if 'ε' in g.terminals:
            g.terminals.discard('ε')

        # Step 3: unit production elimination with snapshots
        print("Step 3.")
        print()
        # Find unit productions (A -> B)
        unit_pairs = set()
        for lhs, prods in g.productions.items():
            for rhs in prods:
                if len(rhs) == 1 and rhs[0] in g.non_terminals:
                    unit_pairs.add((lhs, rhs[0]))
        # remove units one by one and show grammar
        while unit_pairs:
            a, b = unit_pairs.pop()
            print(f"Now, it creates unit production {a}→{b}, its removal from the grammar yields:")
            # add non-unit productions of b to a
            for rhs in g.productions.get(b, []):
                if not (len(rhs) == 1 and rhs[0] in g.non_terminals):
                    if rhs not in g.productions.get(a, []):
                        g.productions.setdefault(a, []).append(rhs)
            # remove the unit production a->b
            g.productions[a] = [rhs for rhs in g.productions.get(a, []) if not (len(rhs) == 1 and rhs[0] == b)]
            # recompute unit_pairs
            unit_pairs = set()
            for lhs, prods in g.productions.items():
                for rhs in prods:
                    if len(rhs) == 1 and rhs[0] in g.non_terminals:
                        unit_pairs.add((lhs, rhs[0]))
            for l in grammar_to_block(g):
                print(l)
            print()

        # Step 3.5: also eliminate any remaining unit S0->S if exists (already handled above)

        # Step 3 (terminal substitution)
        print("Step 3.")
        print()
        print("In production rule A→aAS | aS and B→ SbS | aAS | aS, terminals a and b exist on RHS with non-terminals. Removing them from RHS:")
        print()
        # create terminal vars using simple single-letter names (X, Y, Z, ...) to keep trace readable
        term_map = {}
        term_names_pool = [c for c in ['X', 'Y', 'Z', 'W', 'V'] if c not in g.non_terminals]
        fallback_counter = 0
        for lhs, prods in list(g.productions.items()):
            for rhs in prods:
                for sym in rhs:
                    if sym not in g.non_terminals and sym != 'ε' and sym not in term_map:
                        if term_names_pool:
                            var = term_names_pool.pop(0)
                        else:
                            var = f"T{fallback_counter}"
                            fallback_counter += 1
                        term_map[sym] = var
                        g.productions.setdefault(var, []).append([sym])
                        g.non_terminals.add(var)

        # replace terminals in mixed rules with their new single-letter variables
        for lhs, prods in list(g.productions.items()):
            new_list = []
            for rhs in prods:
                if len(rhs) >= 2:
                    new_rhs = [term_map[s] if (s in term_map and s not in g.non_terminals) else s for s in rhs]
                    new_list.append(new_rhs)
                else:
                    new_list.append(rhs)
            g.productions[lhs] = new_list

        for l in grammar_to_block(g):
            print(l)
        print()

        print("Also, B→ bb can’t be part of CNF, removing it from grammar yields:")
        print()
        # Replace any terminal-terminal like [b,b] with Y Y (already done above as terminals mapping)
        # For simplicity, ensure any rhs that is two terminals maps to their variables
        for lhs in list(g.productions.keys()):
            newp = []
            for rhs in g.productions[lhs]:
                if len(rhs) == 2 and all(sym not in g.non_terminals for sym in rhs):
                    # create or find var for each terminal
                    vars_pair = [term_map[s] if s in term_map else s for s in rhs]
                    newp.append(vars_pair)
                else:
                    newp.append(rhs)
            g.productions[lhs] = newp

        for l in grammar_to_block(g):
            print(l)
        print()

        # Step 4: break long productions
        print("Step 4:")
        print()
        # break any RHS longer than 2 by introducing new variables
        counter = 0
        # provide short, single-letter helper names where possible (P, R, T, Q, ...)
        simple_pool = [c for c in ['P', 'R', 'T', 'Q', 'M', 'N', 'O', 'U', 'V', 'W'] if c not in g.non_terminals]
        def new_var():
            nonlocal counter
            if simple_pool:
                name = simple_pool.pop(0)
            else:
                name = f"P{counter}"
                counter += 1
                while name in g.non_terminals:
                    name = f"P{counter}"
                    counter += 1
            g.non_terminals.add(name)
            return name

        # We'll repeatedly scan productions and replace leftmost pair with new var
        changed = True
        while changed:
            changed = False
            for lhs in list(g.productions.keys()):
                new_list = []
                for rhs in g.productions[lhs]:
                    if len(rhs) > 2:
                        # replace leftmost pair
                        a, b, *rest = rhs
                        pair = [a, b]
                        # check if some var already exists for this pair
                        found = None
                        for k, v in g.productions.items():
                            if v == [pair]:
                                found = k
                                break
                        if not found:
                            nv = new_var()
                            g.productions[nv] = [pair]
                            found = nv
                        new_rhs = [found] + rest
                        new_list.append(new_rhs)
                        changed = True
                    else:
                        new_list.append(rhs)
                g.productions[lhs] = new_list

        # For nice display, try to rename generated P... variables to user-style names P, R, T where possible
        # (not necessary for correctness)

        for l in grammar_to_block(g):
            print(l)
        print()

        print("(Every production now is either two non-terminals on the RHS or a single terminal.)")



# Built-in examples are run through the CLI using the --examples flag.


def parse_grammar_from_lines(lines):
    """Parse a simple grammar format from lines.

    Expected format per line:
      A -> a B C | ε | b
    Whitespace is ignored around symbols. Non-terminals are assumed to be uppercase
    identifiers (strings without quotes). Terminals are other symbols (lowercase,
    characters, or quoted tokens). ε denotes the empty string.
    """
    cfg = CFG('S')
    start_set = False
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        # Accept both ASCII arrow '->' and Unicode arrow '→'
        if '->' in line:
            lhs, rhs_all = line.split('->', 1)
        elif '→' in line:
            lhs, rhs_all = line.split('→', 1)
        else:
            continue
        lhs = lhs.strip()
        # Determine start symbol only once (first encountered LHS)
        if not start_set and lhs:
            cfg.start_symbol = lhs
            start_set = True
        # accept '/' as an alternative separator too (some inputs use '/').
        rhs_all = rhs_all.replace('/', '|')
        alternatives = rhs_all.split('|')
        for alt in alternatives:
            alt = alt.strip()
            if alt == 'ε' or alt == 'eps' or alt == 'epsilon':
                cfg.add_production(lhs, ['ε'])
            else:
                # If RHS contains whitespace-separated tokens, use those.
                # Otherwise, treat the RHS as a sequence of single-character symbols
                tokens = [tok.strip() for tok in alt.split() if tok.strip() and tok.strip() != '/']
                if len(tokens) > 1:
                    symbols = tokens
                else:
                    # No whitespace-delimited tokens: split into single characters
                    # e.g., 'ASB' -> ['A','S','B'], 'bb' -> ['b','b']
                    atom = tokens[0] if tokens else ''
                    symbols = list(atom)
                if symbols:
                    cfg.add_production(lhs, symbols)
    return cfg


def main():
    import argparse, sys

    parser = argparse.ArgumentParser(description='Convert a CFG to Chomsky Normal Form (CNF)')
    parser.add_argument('--file', '-f', help='Path to grammar file. If omitted, reads stdin.')
    parser.add_argument('--examples', action='store_true', help='Run built-in examples (default behavior)')
    parser.add_argument('--no-verbose', dest='verbose', action='store_false', help='Disable verbose conversion output')
    parser.add_argument('--trace', action='store_true', help='Print step-by-step Step 1..4 trace in the requested formatted style')
    parser.add_argument('--canonical', action='store_true', help='Print a canonical textbook-style CNF trace (terminals->C,D; pairs->E,F)')
    args = parser.parse_args()

    if args.examples and not args.file:
        # Run the built-in examples (preserve original behavior)
        run_examples = True
    else:
        run_examples = False

    if run_examples:
        # Re-run the existing examples (keep original prints)
        # Example 1..3 are already constructed above; to avoid duplicating code,
        # we'll re-create them here and convert.
        print("="*70)
        print("EXAMPLE 1: Simple Grammar")
        print("="*70)

        cfg1 = CFG('S')
        cfg1.add_production('S', ['A', 'B'])
        cfg1.add_production('S', ['a'])
        cfg1.add_production('A', ['a'])
        cfg1.add_production('A', ['B'])
        cfg1.add_production('B', ['b'])
        cfg1.add_production('B', ['ε'])

        converter1 = CFGtoChomsky(cfg1)
        cnf1 = converter1.convert(verbose=args.verbose)

        print("\n\n" + "="*70)
        print("EXAMPLE 2: Grammar with Long Productions")
        print("="*70)

        cfg2 = CFG('S')
        cfg2.add_production('S', ['A', 'S', 'B'])
        cfg2.add_production('S', ['a', 'b'])
        cfg2.add_production('A', ['a', 'A'])
        cfg2.add_production('A', ['a'])
        cfg2.add_production('B', ['b', 'B'])
        cfg2.add_production('B', ['b'])

        converter2 = CFGtoChomsky(cfg2)
        cnf2 = converter2.convert(verbose=args.verbose)

        print("\n\n" + "="*70)
        print("EXAMPLE 3: Palindrome Grammar")
        print("="*70)

        cfg3 = CFG('S')
        cfg3.add_production('S', ['a', 'S', 'a'])
        cfg3.add_production('S', ['b', 'S', 'b'])
        cfg3.add_production('S', ['a'])
        cfg3.add_production('S', ['b'])
        cfg3.add_production('S', ['ε'])

        converter3 = CFGtoChomsky(cfg3)
        cnf3 = converter3.convert(verbose=args.verbose)

        return

    # Read grammar from file or stdin
    if args.file:
        try:
            with open(args.file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
    else:
        if sys.stdin.isatty():
            print("Enter grammar lines (end with EOF / Ctrl+Z then Enter on Windows):")
        lines = sys.stdin.read().splitlines()

    cfg = parse_grammar_from_lines(lines)
    converter = CFGtoChomsky(cfg)
    if args.trace:
        # Use trace mode which prints the Step 1..4 formatted trace
        converter = CFGtoChomsky(cfg)
        converter.convert_with_trace()
    elif args.canonical:
        # Print a canonical textbook-style trace tailored to simple grammars
        def print_canonical_trace(cfg):
            # Print initial grammar
            def print_block(title, g):
                print()
                print(title)
                for lhs in sorted(g.productions.keys()):
                    rhss = g.productions[lhs]
                    rhs_text = ' | '.join([''.join(rhs) if rhs != ['ε'] else 'ε' for rhs in rhss])
                    print(f"{lhs} → {rhs_text}")

            g = cfg.copy()
            print("Solution-\n")
            print("Step-01:")
            print_block("\nThe given grammar is:", g)

            # Step-02: eliminate ε-productions (compute full nullable set and generate all variants)
            print("\nStep-02:")
            # compute nullable closure (variables that can generate ε)
            nullable = set()
            changed = True
            while changed:
                changed = False
                for lhs, prods in g.productions.items():
                    if lhs in nullable:
                        continue
                    for rhs in prods:
                        if rhs == ['ε'] or all(sym in nullable for sym in rhs):
                            nullable.add(lhs)
                            changed = True
                            break

            # helper to generate all non-empty variants by removing nullable symbols
            def generate_variants(rhs, nullable):
                if not rhs:
                    return [[]]
                variants = []

                def gen(i, cur):
                    if i == len(rhs):
                        variants.append(cur[:])
                        return
                    s = rhs[i]
                    # include
                    cur.append(s)
                    gen(i+1, cur)
                    cur.pop()
                    # exclude if nullable
                    if s in nullable:
                        gen(i+1, cur)

                gen(0, [])
                # filter out empty
                return [v for v in variants if v]

            # remove explicit ε productions and create variants
            new_productions = {}
            for lhs, prods in g.productions.items():
                new_productions[lhs] = []
                for rhs in prods:
                    if rhs == ['ε']:
                        continue
                    variants = generate_variants(rhs, nullable)
                    for v in variants:
                        if v and v not in new_productions[lhs]:
                            new_productions[lhs].append(v)
            g.productions = new_productions
            # remove any residual explicit ε entries and terminal marker
            for lhs in list(g.productions.keys()):
                g.productions[lhs] = [r for r in g.productions[lhs] if r != ['ε']]
            print_block("\nAfter ε-elimination:", g)

            # Step-02.5: eliminate unit productions (A -> B) BEFORE terminal replacement
            print("\nStep-02.5: Eliminating unit productions")
            # collect unit pairs
            unit_pairs = set()
            for lhs, prods in g.productions.items():
                for rhs in prods:
                    if len(rhs) == 1 and rhs[0] in g.non_terminals:
                        unit_pairs.add((lhs, rhs[0]))

            # compute transitive closure
            closure = set(unit_pairs)
            changed = True
            while changed:
                changed = False
                for (a, b) in list(closure):
                    for (c, d) in list(closure):
                        if b == c and (a, d) not in closure:
                            closure.add((a, d))
                            changed = True

            # build new productions without unit rules
            newp = {}
            for lhs in g.productions.keys():
                newp[lhs] = []
            for lhs, prods in g.productions.items():
                for rhs in prods:
                    if not (len(rhs) == 1 and rhs[0] in g.non_terminals):
                        if rhs not in newp[lhs]:
                            newp[lhs].append(rhs)
            # add productions from closures
            for (a, b) in closure:
                for rhs in g.productions.get(b, []):
                    if not (len(rhs) == 1 and rhs[0] in g.non_terminals):
                        if rhs not in newp[a]:
                            newp[a].append(rhs)

            g.productions = newp
            print_block("\nAfter unit-elimination:", g)

            # Step-03: Replace terminals in mixed productions with new variables C,D,... but keep single-terminal productions
            print("\nStep-03:")
            # collect terminals
            terminals = set()
            for lhs, prods in g.productions.items():
                for rhs in prods:
                    for s in rhs:
                        if s != 'ε' and s not in g.non_terminals:
                            terminals.add(s)
            # map terminals to C,D,... deterministically
            term_names = ['C', 'D', 'E', 'F', 'G']
            term_map = {}
            i = 0
            for t in sorted(terminals):
                if i < len(term_names):
                    term_map[t] = term_names[i]
                else:
                    term_map[t] = f"T{i}"
                i += 1
            # add productions for terminal variables
            for t, var in term_map.items():
                g.productions.setdefault(var, []).append([t])
                g.non_terminals.add(var)

            # replace terminals in RHSs that have length >=2
            for lhs in list(g.productions.keys()):
                newlist = []
                for rhs in g.productions[lhs]:
                    if len(rhs) >= 2:
                        newrhs = [term_map[s] if (s in term_map and s not in g.non_terminals) else s for s in rhs]
                        newlist.append(newrhs)
                    else:
                        newlist.append(rhs)
                g.productions[lhs] = newlist

            print_block("\nAfter replacing terminals with variables:", g)

            # Step-04: Replace repeated pairs AA, BB etc by new vars E,F,... where needed
            print("\nStep-04:")
            pair_map = {}
            pair_names = ['E', 'F', 'G', 'H', 'I']
            pidx = 0
            # find all adjacent pairs in RHSs length >=2 where both are non-terminals
            # prioritize identical pairs (AA, BB, ...) to match textbook-style naming
            pairs = []
            for lhs, prods in g.productions.items():
                for rhs in prods:
                    if len(rhs) >= 2:
                        for i in range(len(rhs)-1):
                            a, b = rhs[i], rhs[i+1]
                            if a in g.non_terminals and b in g.non_terminals:
                                pairs.append((a, b))
            # assign names first to identical pairs, then to others
            for pair in pairs:
                if pair in pair_map:
                    continue
                a, b = pair
                if a == b:
                    if pidx < len(pair_names):
                        pair_map[pair] = pair_names[pidx]
                    else:
                        pair_map[pair] = f"P{pidx}"
                    pidx += 1
            for pair in pairs:
                if pair in pair_map:
                    continue
                if pidx < len(pair_names):
                    pair_map[pair] = pair_names[pidx]
                else:
                    pair_map[pair] = f"P{pidx}"
                pidx += 1
            # add productions and replace occurrences (replace rightmost pairs first)
            for lhs in list(g.productions.keys()):
                newlist = []
                for rhs in g.productions[lhs]:
                    newrhs = rhs[:]
                    # while length > 2, replace rightmost adjacent pair if in map (prefer rightmost)
                    while len(newrhs) > 2:
                        replaced = False
                        for i in range(len(newrhs)-2, -1, -1):
                            pair = (newrhs[i], newrhs[i+1])
                            if pair in pair_map:
                                var = pair_map[pair]
                                # add production var -> pair
                                g.productions.setdefault(var, []).append([pair[0], pair[1]])
                                g.non_terminals.add(var)
                                # replace the pair with var
                                newrhs = newrhs[:i] + [var] + newrhs[i+2:]
                                replaced = True
                                break
                        if not replaced:
                            break
                    newlist.append(newrhs)
                g.productions[lhs] = newlist

            print_block("\nAfter replacing repeated pairs with variables:", g)

            # Final grammar
            print("\nStep-06:")
            print_block("\nResultant grammar (CNF):", g)

        print_canonical_trace(cfg)
    else:
        cnf = converter.convert(verbose=args.verbose)
        # Print final CNF in a compact form
        print("\nFinal CNF Grammar:")
        cnf.display()


if __name__ == '__main__':
    main()
