from expression import predict_from_meta
from expression import preprocess_matrix_for_ct_prediction

########
# usable for unique markers only (different from high/low markers)
########

def celltypes_prediction(cxr_gxc, louvain_arr, gene_names):
	cxg_csr, gene_names = preprocess_matrix_for_ct_prediction(cxr_gxc, gene_names)

	final_ct_prediction_res = predict_from_meta(
		cxg_csr, gene_names, louvain_arr
	)

	return final_ct_prediction_res