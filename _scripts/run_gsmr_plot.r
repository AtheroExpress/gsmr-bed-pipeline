source("vendor/gsmr_plot.r")


args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop("usage: _scripts/run_gsmr_plot.r gsmr_datafile gsmr_filtered plot_dir");
}

gsmr_datafile = args[1]
gsmr_filtered = args[2]
plot_dir = args[3]

gsmr_data = read_gsmr_data(gsmr_datafile) # this will take LONG

df = read.table(gsmr_filtered, header=T)
for (row in 1:nrow(df)) {
    exposure = as.character(df[row, 'Exposure'])
    outcome  = as.character(df[row, 'Outcome'])
    filename = sprintf("%s/%s-%s.pdf", plot_dir, exposure, outcome)
    pdf(file=filename, width=4, height=4)
    plot_gsmr_effect(gsmr_data, exposure, outcome)
    dev.off()
}

