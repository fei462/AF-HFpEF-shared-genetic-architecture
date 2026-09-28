# Optional central environment configuration.
# The released scripts retain the original D:/A/data fallback used for the study,
# but reviewers can relocate the project by defining these environment variables.
Sys.setenv(AFHFPEF_DATA_ROOT = Sys.getenv("AFHFPEF_DATA_ROOT", unset = "D:/A/data"))
# Sys.setenv(LDSC_ROOT = "/path/to/ldsc")
# Sys.setenv(LDSC_PYTHON = "/path/to/python")
