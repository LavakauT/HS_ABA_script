
# ###### available chromosome ######
# peak <- read.delim('/Users/user/Desktop/f1/counts/split/DE_peaks_hclust.bed',
#                    header = F)
# head(peak)
# pass <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chrU', 'chrV',
#           'unplaced-scaffold_010', 'unplaced-scaffold_056', 'unplaced-scaffold_078', 'unplaced-scaffold_086',
#           'unplaced-scaffold_098', 'unplaced-scaffold_131', 'unplaced-scaffold_145', 'unplaced-scaffold_162',
#           'unplaced-scaffold_190', 'unplaced-scaffold_202', 'unplaced-scaffold_221', 'unplaced-scaffold_250',
#           'unplaced-scaffold_257', 'unplaced-scaffold_267', 'unplaced-scaffold_279', 'unplaced-scaffold_281',
#           'unplaced-scaffold_315', 'unplaced-scaffold_349', 'unplaced-scaffold_362', 'unplaced-scaffold_406',
#           'unplaced-scaffold_431', 'unplaced-scaffold_432', 'unplaced-scaffold_433', 'unplaced-scaffold_435',
#           'unplaced-scaffold_436', 'unplaced-scaffold_439', 'unplaced-scaffold_440', 'unplaced-scaffold_441',
#           'unplaced-scaffold_445')
# 
# export <- peak %>% 
#   filter(V1 %in% pass) # 5946 >> 5941
# 
# write.table(export, '/Users/user/Desktop/f1/counts/split/DE_peaks_hclust_available.bed',
#             row.names = F, col.names = F, quote = F, sep = '\t')
# 
# 
# 
# 
# removal <- c('unplaced-scaffold_377', 'unplaced-scaffold_313', 'unplaced-scaffold_038', 'unplaced-scaffold_213',
#              'unplaced-scaffold_405', 'unplaced-scaffold_269', 'unplaced-scaffold_154', 'unplaced-scaffold_160',
#              'unplaced-scaffold_084', 'unplaced-scaffold_016', 'unplaced-scaffold_308', 'unplaced-scaffold_150',
#              'unplaced-scaffold_185', 'unplaced-scaffold_030', 'unplaced-scaffold_089', 'unplaced-scaffold_264',
#              'unplaced-scaffold_375', 'unplaced-scaffold_359', 'unplaced-scaffold_115', 'unplaced-scaffold_385',
#              'unplaced-scaffold_268', 'unplaced-scaffold_383', 'unplaced-scaffold_307', 'unplaced-scaffold_410',
#              'unplaced-scaffold_156', 'unplaced-scaffold_059', 'unplaced-scaffold_381', 'unplaced-scaffold_153',
#              'unplaced-scaffold_258', 'unplaced-scaffold_358', 'unplaced-scaffold_425', 'unplaced-scaffold_230',
#              'unplaced-scaffold_044', 'unplaced-scaffold_014', 'unplaced-scaffold_079', 'unplaced-scaffold_209',
#              'unplaced-scaffold_168', 'unplaced-scaffold_380', 'unplaced-scaffold_397', 'unplaced-scaffold_320',
#              'unplaced-scaffold_335', 'unplaced-scaffold_283', 'unplaced-scaffold_126', 'unplaced-scaffold_282',
#              'unplaced-scaffold_379', 'unplaced-scaffold_338', 'unplaced-scaffold_233', 'unplaced-scaffold_139',
#              'unplaced-scaffold_325', 'unplaced-scaffold_271', 'unplaced-scaffold_403', 'unplaced-scaffold_028',
#              'unplaced-scaffold_161', 'unplaced-scaffold_108', 'unplaced-scaffold_295', 'unplaced-scaffold_149',
#              'unplaced-scaffold_085', 'unplaced-scaffold_340', 'unplaced-scaffold_092', 'unplaced-scaffold_183',
#              'unplaced-scaffold_395', 'unplaced-scaffold_081', 'unplaced-scaffold_029', 'unplaced-scaffold_065',
#              'unplaced-scaffold_345', 'unplaced-scaffold_104', 'unplaced-scaffold_072', 'unplaced-scaffold_187',
#              'unplaced-scaffold_054', 'unplaced-scaffold_260', 'unplaced-scaffold_046', 'unplaced-scaffold_199',
#              'unplaced-scaffold_322', 'unplaced-scaffold_404', 'unplaced-scaffold_101', 'unplaced-scaffold_091',
#              'unplaced-scaffold_095', 'unplaced-scaffold_058', 'unplaced-scaffold_222', 'unplaced-scaffold_111',
#              'unplaced-scaffold_076', 'unplaced-scaffold_170', 'unplaced-scaffold_147', 'unplaced-scaffold_164',
#              'unplaced-scaffold_418', 'unplaced-scaffold_179', 'unplaced-scaffold_060', 'unplaced-scaffold_148',
#              'unplaced-scaffold_417', 'unplaced-scaffold_011', 'unplaced-scaffold_212', 'unplaced-scaffold_354',
#              'unplaced-scaffold_146', 'unplaced-scaffold_408', 'unplaced-scaffold_067', 'unplaced-scaffold_171',
#              'unplaced-scaffold_259', 'unplaced-scaffold_073', 'unplaced-scaffold_270', 'unplaced-scaffold_394',
#              'unplaced-scaffold_158', 'unplaced-scaffold_239', 'unplaced-scaffold_292', 'unplaced-scaffold_372',
#              'unplaced-scaffold_018', 'unplaced-scaffold_364', 'unplaced-scaffold_032', 'unplaced-scaffold_035',
#              'unplaced-scaffold_195', 'unplaced-scaffold_365', 'unplaced-scaffold_129', 'unplaced-scaffold_192',
#              'unplaced-scaffold_424', 'unplaced-scaffold_275', 'unplaced-scaffold_414', 'unplaced-scaffold_013',
#              'unplaced-scaffold_353', 'unplaced-scaffold_045', 'unplaced-scaffold_303', 'unplaced-scaffold_134',
#              'unplaced-scaffold_420', 'unplaced-scaffold_291', 'unplaced-scaffold_336', 'unplaced-scaffold_042',
#              'unplaced-scaffold_430', 'unplaced-scaffold_075', 'unplaced-scaffold_319', 'unplaced-scaffold_140',
#              'unplaced-scaffold_261', 'unplaced-scaffold_142', 'unplaced-scaffold_074', 'unplaced-scaffold_272',
#              'unplaced-scaffold_419', 'unplaced-scaffold_262', 'unplaced-scaffold_280', 'unplaced-scaffold_242',
#              'unplaced-scaffold_248', 'unplaced-scaffold_407', 'unplaced-scaffold_173', 'unplaced-scaffold_197', 
#              'unplaced-scaffold_031', 'unplaced-scaffold_226', 'unplaced-scaffold_122', 'unplaced-scaffold_055',
#              'unplaced-scaffold_297', 'unplaced-scaffold_331', 'unplaced-scaffold_136', 'unplaced-scaffold_097',
#              'unplaced-scaffold_124', 'unplaced-scaffold_121', 'unplaced-scaffold_266', 'unplaced-scaffold_087',
#              'unplaced-scaffold_423', 'unplaced-scaffold_373', 'unplaced-scaffold_062', 'unplaced-scaffold_413',
#              'unplaced-scaffold_041', 'unplaced-scaffold_216', 'unplaced-scaffold_157', 'unplaced-scaffold_309',
#              'unplaced-scaffold_330', 'unplaced-scaffold_159', 'unplaced-scaffold_241', 'unplaced-scaffold_400',
#              'unplaced-scaffold_177', 'unplaced-scaffold_246', 'unplaced-scaffold_205', 'unplaced-scaffold_391',
#              'unplaced-scaffold_019', 'unplaced-scaffold_247', 'unplaced-scaffold_387', 'unplaced-scaffold_351',
#              'unplaced-scaffold_399', 'unplaced-scaffold_274', 'unplaced-scaffold_021', 'unplaced-scaffold_052',
#              'unplaced-scaffold_426', 'unplaced-scaffold_256', 'unplaced-scaffold_228', 'unplaced-scaffold_070', 
#              'unplaced-scaffold_316', 'unplaced-scaffold_285', 'unplaced-scaffold_143', 'unplaced-scaffold_276',
#              'unplaced-scaffold_255', 'unplaced-scaffold_341', 'unplaced-scaffold_061', 'unplaced-scaffold_204',
#              'unplaced-scaffold_165', 'unplaced-scaffold_198', 'unplaced-scaffold_402', 'unplaced-scaffold_163', 
#              'unplaced-scaffold_286', 'unplaced-scaffold_231', 'unplaced-scaffold_376', 'unplaced-scaffold_144',
#              'unplaced-scaffold_304', 'unplaced-scaffold_043', 'unplaced-scaffold_117', 'unplaced-scaffold_284',
#              'unplaced-scaffold_339', 'unplaced-scaffold_012', 'unplaced-scaffold_099', 'unplaced-scaffold_050',
#              'unplaced-scaffold_327', 'unplaced-scaffold_224', 'unplaced-scaffold_203', 'unplaced-scaffold_298',
#              'unplaced-scaffold_232', 'unplaced-scaffold_263', 'unplaced-scaffold_273', 'unplaced-scaffold_127',
#              'unplaced-scaffold_235', 'unplaced-scaffold_053', 'unplaced-scaffold_109', 'unplaced-scaffold_389',
#              'unplaced-scaffold_278', 'unplaced-scaffold_210', 'unplaced-scaffold_392', 'unplaced-scaffold_444',
#              'unplaced-scaffold_314', 'unplaced-scaffold_094', 'unplaced-scaffold_393', 'unplaced-scaffold_103',
#              'unplaced-scaffold_141', 'unplaced-scaffold_390', 'unplaced-scaffold_305', 'unplaced-scaffold_034',
#              'unplaced-scaffold_066', 'unplaced-scaffold_411', 'unplaced-scaffold_288', 'unplaced-scaffold_384')
# 
# 
# peak <- read.delim('/Users/user/Desktop/f1/counts/split/DE_peaks_hclust.bed',
#                    header = F)
# export <- peak %>% 
#   filter(!(V1 %in% removal)) # 5946 >> 5941
# 
# write.table(export, '/Users/user/Desktop/f1/counts/split/DE_peaks_hclust_available.bed',
#             row.names = F, col.names = F, quote = F, sep = '\t')
# ##############################
