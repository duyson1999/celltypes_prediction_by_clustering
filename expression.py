import json
import numpy as np
from sklearn.utils import sparsefuncs
from constants import DEFAULT_LIBRARY_SIZE_FACTOR
from example_markers import CELLTYPES_PREDICTION as T_HELPER2_ALLERGY_MARKERS


def _init_info():
    return {
        'expressed_coverage': 0,
        'entropy': 0,
        'mean_expression': 0,
        'quantile_50': 0
    }


def _convert_boolean_exp(expression_arr, cutoff):
    bin_mtx = np.zeros_like(expression_arr)
    bin_mtx[expression_arr > cutoff] = 1
    return bin_mtx


def _convert_library_size(csr_cxg, size_factor=DEFAULT_LIBRARY_SIZE_FACTOR):
    res = []
    for i in range(50):
        res.append(
            np.sum(
                np.expm1(
                    csr_cxg.data[csr_cxg.indptr[i] : csr_cxg.indptr[i+1]]
                )
            )
        )
    res = [np.round(i) for i in res]
    ori_size_factor = np.unique(res)
    upper_bound = ori_size_factor[0] + ori_size_factor[0] * 0.1
    lower_bount = ori_size_factor[0] - ori_size_factor[0] * 0.1

    if np.all(ori_size_factor > lower_bount) and np.all(ori_size_factor < upper_bound):
        if size_factor < upper_bound and size_factor > lower_bount:
            return csr_cxg

        csr_cxg.data = np.expm1(csr_cxg.data)
        csr_cxg.data = csr_cxg.data * (size_factor / ori_size_factor[0])
        csr_cxg.data = np.log1p(csr_cxg.data)
        return csr_cxg

    print('cannot find the library size:', ori_size_factor)
    return None


def _calculate_stats(cxg_csr, cutoff):
    expression_arr = cxg_csr.T.toarray().reshape(-1)
    # print (expression_arr)
    bin_mtx = _convert_boolean_exp(expression_arr, cutoff)
    res = _init_info()
    res['expressed_coverage'] = np.sum(bin_mtx > 0) / len(bin_mtx)
    res['entropy'] = - (res['expressed_coverage'] * np.log2(res['expressed_coverage']))
    res['mean_expression'] = np.sum(expression_arr[expression_arr > cutoff]) / len(expression_arr)
    res['quantile_50'] = np.quantile(expression_arr, 0.5)
    return res


def _get_basic_gene_exp_info(cxg_csr, gene_idx):
    all_cells_exp = cxg_csr[:, gene_idx].data
    if not len(all_cells_exp):
        return 0, 0
    # cutoff = 0.1 * all_cells_exp.max()
    #
    tmp_exp = np.quantile(all_cells_exp, 0.95)
    tmp_val = np.sort(all_cells_exp)[all_cells_exp < tmp_exp][-1]
    cutoff = 0.4 * tmp_val
    #
    all_cells_mean_expression = np.sum(all_cells_exp[all_cells_exp > cutoff]) / cxg_csr.shape[0]
    return cutoff, all_cells_mean_expression


def _cluster_express_gene(
        cxg_csr,
        gene_idx,
        cluster_bool_arr,
        cutoff,
        all_cells_mean_expression
    ):
    tmp_res = _calculate_stats(cxg_csr[np.nonzero(cluster_bool_arr)[0], :][:, gene_idx], cutoff)
    if tmp_res['expressed_coverage'] >= 0.5:
        return True
    else:
        return False
    # if tmp_res['expressed_coverage'] < 0.1:
    #     return False
    # elif tmp_res['expressed_coverage'] >= 0.7:
    #     return True
    # elif tmp_res['expressed_coverage'] < 0.7 \
    #     and tmp_res['mean_expression'] < all_cells_mean_expression:
    #     return False
    # else:
    #     return True


def _precalculate_genes_info(cxg_csr, gene_names_annotation, all_genes):
    all_genes = sorted(list(set(all_genes)))
    genes_info = {
        i: {
            'cutoff': 0,
            'all_cells_expressed_coverage': 0,
            'all_cells_mean_expression': 0
        } for i in all_genes
    }
    for gene in all_genes:
        print (gene)
        gene_idx = np.nonzero(gene_names_annotation == gene)[0][0]
        x, y = _get_basic_gene_exp_info(
            cxg_csr,
            gene_idx
        )
        print ('cutoff', x)
        print ('all_cells_mean_expression', y)
        genes_info[gene]['cutoff'] = x
        genes_info[gene]['all_cells_mean_expression'] = y
    return genes_info


def _precalculate_expressed_clusters(
        cxg_csr,
        gene_names_annotation,
        louvain_arr,
        genes_info,
        all_genes
    ):
    genes_cluster_bool = {
        i: np.zeros(len(louvain_arr), dtype=np.bool_)
        for i in all_genes
    }
    for gene in all_genes:
        x = genes_info[gene]['cutoff']
        y = genes_info[gene]['all_cells_mean_expression']
        gene_idx = np.nonzero(gene_names_annotation == gene)[0][0]
        for cluster_name in np.unique(louvain_arr):
            cluster_bool_arr = louvain_arr == cluster_name
            is_expressed = _cluster_express_gene(
                cxg_csr,
                gene_idx,
                cluster_bool_arr,
                x,
                y,
            )
            genes_cluster_bool[gene][cluster_bool_arr] = is_expressed
    return genes_cluster_bool


def _precalculate_prediction(
        cxg_csr,
        gene_names_annotation,
        louvain_arr,
        context,
        markers_list=T_HELPER2_ALLERGY_MARKERS
    ):
    all_genes = []
    for i in markers_list:
        for j in markers_list[i]:
            all_genes.extend(j[0])
            all_genes.extend(j[1])
    all_genes = sorted(list(set(all_genes)))

    genes_info = _precalculate_genes_info(cxg_csr, gene_names_annotation, all_genes)

    if len(context):
        cxg_csr = cxg_csr[context, :]
        louvain_arr = louvain_arr[context]
    genes_cluster_bool = _precalculate_expressed_clusters(
        cxg_csr,
        gene_names_annotation,
        louvain_arr,
        genes_info,
        all_genes
    )

    tmp_save = genes_cluster_bool.copy()
    for i in tmp_save:
        tmp_save[i] = np.unique(louvain_arr[tmp_save[i]]).tolist()
    with open('/data/sonvo/atlas/expression_log.json', 'w') as f:
        json.dump(tmp_save, f)

    return genes_cluster_bool


def _define_celltypes(
        n_cells,
        gene_lists,
        genes_cluster_bool
    ):
    positive_tmp_bool_arr = np.zeros(shape=n_cells, dtype=np.bool_)
    if len(gene_lists[0]) == 0:
        positive_tmp_bool_arr[:] = True
    for gene_name in gene_lists[0]:
        is_expressed_arr = genes_cluster_bool[gene_name]
        positive_tmp_bool_arr[is_expressed_arr] = True

    negative_tmp_bool_arr = np.ones(shape=n_cells, dtype=np.bool_)
    for gene_name in gene_lists[1]:
        is_expressed_arr = genes_cluster_bool[gene_name]
        negative_tmp_bool_arr[is_expressed_arr] = False

    return np.logical_and(positive_tmp_bool_arr, negative_tmp_bool_arr)


def _normalize_total(csr_matrix, counts):
	csr_matrix = csr_matrix.astype(np.float32)
	counts += counts == 0
	counts = counts / 10000
	sparsefuncs.inplace_row_scale(csr_matrix, 1 / counts)
	return csr_matrix


def _log_normalize(csr_matrix):
	counts_per_cell = csr_matrix.sum(1)
	counts_per_cell = np.ravel(counts_per_cell)

	csr_matrix = _normalize_total(csr_matrix, counts_per_cell)
	csr_matrix.data = np.log1p(csr_matrix.data)

	return csr_matrix


def preprocess_matrix_for_ct_prediction(cxr_gxc, fe):
	mtx = cxr_gxc.astype(np.float32)
	raw = mtx.data

	ADT_genes = [i for i in range(len(fe)) if fe[i].startswith('ADT-')]
	if len(ADT_genes) != 0:
		genes_bool = np.ones(len(fe), dtype=np.bool_)
		genes_bool[ADT_genes] = False
		mtx = mtx[np.nonzero(genes_bool)[0], :]
		fe = fe[genes_bool]

	mtx = mtx.T.tocsr()

	if (np.sum(raw > 30) > 0 and not np.all(raw.astype('int') == raw)) or np.all(raw.astype('int') == raw):
		mtx = _log_normalize(mtx)
	else:
		mtx = _convert_library_size(mtx)

	return mtx, fe


def predict_from_meta(cxg_csr, fe, louvain_arr, context, level_3_markers=T_HELPER2_ALLERGY_MARKERS):
    n_cells = cxg_csr.shape[0]
    if len(context):
        n_cells = len(context)

    ct_prediction_res = {
        i: np.array([''], dtype='object')[np.zeros(shape=n_cells, dtype=np.uint8)] for i in level_3_markers
    }
    genes_cluster_bool = _precalculate_prediction(cxg_csr, fe, louvain_arr, context, level_3_markers)

    for celltype in level_3_markers.keys():
        all_level_gene_lists = level_3_markers[celltype]
        res = None
        for gene_list in all_level_gene_lists:
            tmp_res = _define_celltypes(
                n_cells,
                gene_list,
                genes_cluster_bool
            )
            if res is None:
                res = tmp_res
                continue
            res = np.logical_and(res, tmp_res)
        ct_prediction_res[celltype][res] = celltype

    final_ct_prediction_res = np.array([''], dtype='object')[
        np.zeros(shape=n_cells, dtype=np.uint8)
    ]
    for i in ct_prediction_res.keys():
        ct_prediction_res[i][ct_prediction_res[i] != ''] = ' --> ' + ct_prediction_res[i][ct_prediction_res[i] != '']

    for i in ct_prediction_res.keys():
        final_ct_prediction_res += ct_prediction_res[i]

    return final_ct_prediction_res

