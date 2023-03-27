
# # # # # # # # # # # # # # # #
# IMPORTS
# # # # # # # # # # # # # # # #

import abagen
from neuromaps.nulls import alexander_bloch
from neuromaps.stats import compare_images
import nibabel as nib
import pandas as pd
from statsmodels.stats.multitest import multipletests

# # # # # # # # # # # # # # # #
# FILES
# # # # # # # # # # # # # # # #

EXPRESSION = 'abagen_expression_dkt.csv'
COMPONENT_TABLE = 'ptcs_prepped.csv'
TEST_COMPONENTS = [
    'Orbitofrontal',
    'LateralFrontal',
    'MedialTemporal',
    'Precuneus',
    'LeftParietalTemporal',
    'Sensorimotor',
    'Occipital',
    'RightParietalTemporal'
    ]
GENE_LIST = 'wightman_genes.txt'


# # # # # # # # # # # # # # # #
# PARAMETERS
# # # # # # # # # # # # # # # #

PROGRESS = 10
N_PERM = 5000
P_CORRECT = ['fdr_bh']

# # # # # # # # # # # # # # # #
# OVERTURE
# # # # # # # # # # # # # # # #

print()
print('NEUROMAPS FOR TAU ANALYSIS')
print('--------------------------')
print(f"Expression data: {EXPRESSION}")
print(f"Component table: {COMPONENT_TABLE}")
print(f"Components tested: {TEST_COMPONENTS}")
print(f"Permutations: {N_PERM}")
print(f"p-value corrections: {P_CORRECT}")

# # # # # # # # # # # # # # # #
# READ DATA
# # # # # # # # # # # # # # # #

print()
print('.......')
print('*** Reading DKT Atlas ***')
dkt = abagen.fetch_desikan_killiany(surface=True)
dkt_l_path, dkt_r_path = dkt['image']
dkt_info_path = dkt['info']

dkt_l = nib.load(dkt_l_path)
dkt_r =  nib.load(dkt_r_path)
parcellation = (dkt_l, dkt_r)

dkt_labels = pd.read_csv(dkt_info_path)
dkt_labels = dkt_labels[dkt_labels['structure'] == 'cortex']
print('SUCCESS')
print('.......')

# read expression data
# could be loaded directly from abagen,
# but for now loading explicitly
print()
print('.......')
print('*** Reading saved expression data ***')
expression = pd.read_csv(EXPRESSION)
genes_only = expression.drop('label', axis=1)
print('SUCCESS')
print('.......')

# read the NMF components,
# already with some preprocessing applied for this purpose
print()
print('.......')
print('*** Reading prepped component data ***')
components = pd.read_csv(COMPONENT_TABLE)
print('SUCCESS')
print('.......')

print()
print('.......')
print('*** Reporting ordering of regions ***')

order_components = components[['label', 'hemisphere']]
order_genes = (dkt_labels.
    set_index('id').
    reindex(expression['label']).
    reset_index(drop=True))[['label', 'hemisphere']]

pd.set_option('display.max_rows', 100)

print()
print("COMPONENTS:")
print(order_components)

print()
print("GENES:")
print(order_genes)

print()
print("ALL EQUAL?")
print((order_components == order_genes).all(axis=None))

print('.......')

# # # # # # # # # # # # # # # #
# READ GENES
# # # # # # # # # # # # # # # #

print()
print('.......')
print('*** Reading gene comparison list ***')
with open(GENE_LIST, 'r') as f:
    genes = [g.strip() for g in f.readlines()]

print()
print('Reporting genes:')
for g in genes:
    print(g)

print()
print('SUCCESS')
print('.......')

# # # # # # # # # # # # # # # #
# MAIN
# # # # # # # # # # # # # # # #

print()
print('!!! BEGINNING MAIN SCRIPT !!!')

for n, component in enumerate(TEST_COMPONENTS):
    print()
    print(f'Component = {component}')

    rows = []

    # this is originally how I set rotations
    # but John brought up issue that using same
    # random seed for all comparisons might have issues

    # now, I instead make a new null for each comparison
    # this takes a bit longer but not terribly so

    # alternative would be to set `seed=n` here and
    # have one null for each component.  I tried
    # all these approaches and found similar significance
    # results - only changes in some associations
    # near the alpha level (which might just be due
    # to differences in seed selection rather than problems
    # with the correlated nulls)

    # rotated = alexander_bloch(components[component],
    #                           atlas='fsaverage',
    #                           density='10k',
    #                           n_perm=N_PERM,
    #                           seed=42,
    #                           parcellation=parcellation)

    for i, gene in enumerate(genes):
        if i % PROGRESS == 0 and i != 0:
            print(f"Gene {i} ({gene})...")

        if gene not in expression.columns:
            print(f"Skipping gene={gene}, not found in expression data...")
            continue

        rotated = alexander_bloch(components[component],
                                  atlas='fsaverage',
                                  density='10k',
                                  n_perm=N_PERM,
                                  seed=8 * n + i,
                                  parcellation=parcellation)

        corr, pval = compare_images(components[component],
                                    expression[gene],
                                    nulls=rotated)

        rows.append({'gene': gene, 'corr': corr, 'pval': pval})

    results = pd.DataFrame(rows)
    for method in P_CORRECT:
        results[method] = multipletests(results['pval'], method=method)[1]

    savename = f"genecorrelation_{component}.csv"
    results.to_csv(savename, index=False)
