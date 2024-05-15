BLANK = []

immune_cells = [
	[['PTPRC'], BLANK]
]
immune_cells__T_cells = immune_cells + [
	[['CD3D', 'CD3G', 'CD3E'], BLANK]
]

_CD4_positive = ['CD4']
_CD8_positive = ['CD8A', 'CD8B']
_gamma_delta_positive = ['TRDC', 'TRDV1', 'TRDV2', 'TRGC1', 'TRGC2', 'TRGV9']
immune_cells__T_cells__CD8_T_cells = immune_cells__T_cells + [
	[_CD8_positive, _gamma_delta_positive + _CD4_positive]
]
immune_cells__T_cells__gamma_delta_T_cells = immune_cells__T_cells + [
	[_gamma_delta_positive, _CD8_positive + _CD4_positive]
]
immune_cells__T_cells__CD4_T_cells = immune_cells__T_cells + [
	[BLANK, _CD8_positive + _gamma_delta_positive]
]

_T_helper_1_positive = ['IFNG']
_T_helper_17_positive = ['IL17A', 'IL17F']
_T_helper_2_positive = ['IL4', 'IL5', 'IL9', 'IL13', 'GATA3-AS1', 'IL17RB', 'PTGDR2', 'HPGDS', 'IL1RL1']
_T_helper_22or2_positive = ['IL22']
_T_reg_positive = ['FOXP3']
_T_folicular_positive = ['CXCR5']
_T_naive = ['CCR7', 'SELL']
immune_cells__T_cells__CD4_T_cells__naive_T_cells = immune_cells__T_cells__CD4_T_cells + [
	[
		_T_naive,
		_T_helper_1_positive + _T_helper_2_positive + _T_helper_17_positive + _T_helper_22or2_positive + _T_reg_positive + _T_folicular_positive
	]
]
# Helper T cells during transitional state could express CCR7 and SELL
HIDDEN_helper_T_cells = immune_cells__T_cells__CD4_T_cells + [
	[
		_T_helper_1_positive + _T_helper_2_positive + _T_helper_17_positive + _T_helper_22or2_positive + _T_reg_positive + _T_folicular_positive,
		BLANK
	]
]

immune_cells__T_cells__CD4_T_cells__T_helper_1 = HIDDEN_helper_T_cells + [
	[
		_T_helper_1_positive,
		_T_helper_2_positive + _T_helper_17_positive + _T_helper_22or2_positive + _T_reg_positive + _T_folicular_positive
	]
]
immune_cells__T_cells__CD4_T_cells__T_helper_17 = HIDDEN_helper_T_cells + [
	[
		_T_helper_17_positive,
		_T_helper_1_positive + _T_helper_2_positive + _T_reg_positive + _T_folicular_positive
	]
]
immune_cells__T_cells__CD4_T_cells__T_helper_2 = HIDDEN_helper_T_cells + [
	[
		_T_helper_2_positive,
		_T_helper_1_positive + _T_helper_17_positive + _T_reg_positive + _T_folicular_positive
	]
]
immune_cells__T_cells__CD4_T_cells__T_helper_22or2 = HIDDEN_helper_T_cells + [
	[
		_T_helper_22or2_positive,
		_T_helper_1_positive + _T_helper_2_positive + _T_helper_17_positive + _T_reg_positive + _T_folicular_positive
	]
]
immune_cells__T_cells__CD4_T_cells__T_reg = HIDDEN_helper_T_cells + [
	[
		_T_reg_positive,
		_T_helper_1_positive + _T_helper_2_positive + _T_helper_22or2_positive + _T_helper_17_positive + _T_folicular_positive
	]
]
immune_cells__T_cells__CD4_T_cells__T_folicular_helper = HIDDEN_helper_T_cells + [
	[
		_T_folicular_positive,
		_T_helper_1_positive + _T_helper_2_positive + _T_helper_22or2_positive + _T_reg_positive + _T_helper_17_positive
	]
]

CELLTYPES_PREDICTION = {
	'immune_cells': immune_cells,
	'T_cells': immune_cells__T_cells,
	'CD8_T_cells': immune_cells__T_cells__CD8_T_cells,
	'gamma_delta_T_cells': immune_cells__T_cells__gamma_delta_T_cells,
	'CD4_T_cells': immune_cells__T_cells__CD4_T_cells,
	'naive_T_cells': immune_cells__T_cells__CD4_T_cells__naive_T_cells,
	'T_helper1': immune_cells__T_cells__CD4_T_cells__T_helper_1,
	'T_helper2': immune_cells__T_cells__CD4_T_cells__T_helper_2,
	'T_helper22/2': immune_cells__T_cells__CD4_T_cells__T_helper_22or2,
	'T_helper17': immune_cells__T_cells__CD4_T_cells__T_helper_17,
	'regulatory_T_cells': immune_cells__T_cells__CD4_T_cells__T_reg,
	'T_folicular_helper': immune_cells__T_cells__CD4_T_cells__T_folicular_helper,
}
