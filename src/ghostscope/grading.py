"""Evidence tiers are conditional routes, not posterior probabilities of ghost ancestry."""
ALL=['Ghost pulse','Direct ingroup pulse','Rate-heterogeneous pure tree','Unexamined histories']

def grade(result):
    flags=list(result.get('flags',[]));s=result.get('spine',{});four=result.get('four_leaf',{});six=result.get('six_leaf',{})
    multi=result.get('multi_haplotype',{});depth=result.get('distance_calibration',{})
    tier='D';title='Model-tied / ambiguous';candidates=list(ALL)
    if s.get('type_III_possible'):
        flags.append('TYPE_III_PARALLEL_DONOR_WARNING')
        next_step='Sample lineages that resolve donor stem/crown placement; compare calendar-calibrated competing models.'
    elif six.get('passed') or (four.get('passed') and result.get('n_ingroup')==4 and not result.get('omitted_sampled_taxa')):
        tier='A';title='Conditional multi-taxon topological support'
        candidates=['Ghost pulse within the certified family','Unexamined histories (multiple pulses, mixtures, ancestral structure, gene-tree bias)']
        next_step='Check backbone, candidate coverage and gene-tree error with independent inference and held-out null simulations.'
        if four.get('passed'):flags.append('UNCALIBRATED_GTEE_TOLERANCE')
    elif multi.get('tier_b_eligible'):
        tier='B';title='Conditional multi-haplotype sufficient-zone support'
        candidates=['Ghost pulse within canonical family','Unexamined histories']
        next_step='Verify common per-locus rate, exchangeability and canonical topology; validate distance measurement error.'
    elif depth.get('quasi_clock_support'):
        tier='C';title='Conditional quasi-clock consistency signal'
        flags.append('QUASI_CLOCK_SCREEN_NOT_CLOCK_PROOF')
        next_step='Fit clock and nonclock competing histories; check ancestral structure and locus-rate mixtures.'
    else:
        next_step='Add a recipient sister and a pure outgroup for rooted four-leaf contrasts; sample outer ladder clades for six-leaf contrasts.'
        if result.get('n_ingroup',0)>=4:
            next_step='Increase independent resolved loci and validate gene-tree error; inspect the weakest retained direct-pulse candidate and backbone assumptions.'
        if 'INSUFFICIENT_TERMINAL_SPINE' in flags:
            next_step='This recipient lacks the certified terminal ladder. Supply a justified alternative backbone or fit competing full networks.'
        if depth.get('sign_flip'):
            tier='E';title='Depth signal unsupported after path calibration';flags.append('SIGN_FLIP')
        if result.get('n_trees',1)==0 and not result.get('multi_haplotype'):
            tier='E';title='No usable topology evidence';next_step='Use --nj-trees or supply independently inferred gene trees.'
    if result.get('triplet',{}).get('open_wedge'):flags.append('SCOPE_ROUTING: OPEN_TOPOLOGY_WEDGE')
    if result.get('six_moment',{}).get('compatible'):flags.append('SCOPE_ROUTING: TWO_RATE_PURE_TREE_MOMENT_COMPATIBLE')
    if multi.get('code'):flags.append(('SCOPE_ROUTING: ' if multi.get('status')=='SCOPE_ROUTING' else '')+multi['code'])
    if result.get('rooting')=='unrooted' and result.get('n_ingroup',0)==4:flags.append('UNROOTED_FOUR_LEAF_FLIP_EQUIVALENCE')
    return {'tier':tier,'title':title,'compatible_candidate_histories':candidates,'flags':sorted(set(flags)),
            'next_step':next_step,'scope':'Evidence tier is conditional on its stated candidate family and observational assumptions; no tier proves the unique biological history.'}
