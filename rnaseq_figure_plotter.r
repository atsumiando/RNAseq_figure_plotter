#Please type this code 'install.packages("library")' to install three libraies 'ggplot2', 'pheatmap', 'tidyverse', and 'argparse' to run rnaseq figure plotter software.

#measure time
start.time <- Sys.time()



#parameter settings
#argparse library
library("argparse")

#parameter setting
parser <- ArgumentParser(usage='Rscript rnaseq_figure_plotter.r -i input_file -t bar -o output_file -g gene_list_file  ... -c 5 -s 6', add_help=TRUE)

#required parameter
parser$add_argument('-i','--input', help = 'input file name', required = TRUE)
parser$add_argument('-t','--type', help = 'choose plot types (bar, box, dot_color, dot_shape, density, density_fill, heatmap, histogram, line, scatter, or violin)', required = TRUE)

#general optional parameter
parser$add_argument('-o', '--output', help = 'default output; output file name',default = "output")
parser$add_argument('-g', '--gene',help = 'default None; file name of specific gene ID list; generate "output"_gene_selection.txt file', default = " ")
parser$add_argument('-r', '--remove_col',help = 'default None; remove specific columns (samples) from input file. Split column name by space. Example; sample1 sample2 sample3', default = " ", nargs= '+')
parser$add_argument('-l', '--log',help = 'default None(0); calculate log value (log2; 2, log10; 10, loge; e)',default = "0")
parser$add_argument('-lgn', '--log_number',help = 'default 0.000000001; add number to avoid -inf for log value',default = 0.000000001, type = "double")
parser$add_argument('-zs', '--zscore',help = 'default off; apply Z-score transformation in gene (on or off). --log function should be 0 to apply --zscore function.',default = "off")
parser$add_argument('-x', '--xaxis',help = 'default samples; choose x-axis (gene or sample)',default = "sample")
parser$add_argument('-z', '--zaxis',help = 'default gene; choose fill, color, or shape (gene or sample)',default = "gene")
parser$add_argument('-c', '--color',help = 'default 1; choose color type (0-10)', default = 1, type = "integer")
parser$add_argument('-cst', '--custom_color',help = 'default None; customize color scales. Split colors by space. Example; red white blue green yellow', default = "0", nargs= '+')
parser$add_argument('-ls', '--letter_size',help = 'default 8 10; type text and title size of legend and axis, respectively. Split two number by space. Example; 20 24', default = c(8,10),type = "double", nargs= 2)
parser$add_argument('-f', '--figure_save_format',help = 'default pdf; choose format of figures (eps, ps, tex (pictex), pdf, jpeg, tiff, png, bmp, svg)', default = 'pdf')

#optional parameter for individual plot types
parser$add_argument('-s', '--style',help = 'default 4; choose backgroud of figures (1-7). This function works for any plots excepts heatmap.', default = 4, type = "integer")
parser$add_argument('-lim', '--limit',help = 'default None; apply individual scale of “data”. This function works for any plots except heatmap. Split two numbers(e.g. limit 0 to 200 -> type 0 200) by space. Negative number required double quotation marks such as “negative number”.  Example; 0 100/“-1” 3', action = 'store', type = "double", nargs= 2, default = c(0,0))
parser$add_argument('-a', '--axis_change',help = 'default off; flip axis in figures (on or off). This function works for any plots except heatmap.', default = "off")
parser$add_argument('-lp', '--legend_position',help = 'default right; choose legend position of figures (none, left, right, bottom, top, or two-element numeric vector). This function works for any plots except heatmap and scatter.', default = "right")
parser$add_argument('-gp', '--geom_position',help = 'default 1; choose visualize types (geom position) from 1-4 in bar, density, density_fill, and histogram', default = 1, type = "integer")
parser$add_argument('-cs', '--cluster_select',help = 'default on on; apply column and row cluster function for heatmap (on or off). Column is first and row is second, split two factor(on or off) by space. Example; on off',default = c('on','on'), nargs= 2)
parser$add_argument('-ss', '--scatter_select',help = 'default None; type column of two samples for comparison in dot plot. Split samples by space. Example; sample1 sample2', nargs= 2)
parser$add_argument('-p', '--plot_size',help = 'default 7 7; type width and height of figure. Split two number by space. This function works for any plots excepts heatmap. Example; 10 12 ', default = c(7,7), type = "double", nargs= 2)

#parameter combine
args <- parser$parse_args()



#parameter name setting
input_file_name <-args$input
output_file_name <-args$output
list_file_name <- args$gene
xaxis <- args$xaxis
zaxis <- args$zaxis
color <- args$color
save_format <- args$figure_save_format



#libraries
library("ggplot2")
library("pheatmap")
library("tidyverse")



#version
version <- "0.0.3"
paste("You are using rnaseq_figure_plotter(R) version; ", version)



#final file name generation
name_setting <-  c(output_file_name, save_format)
collapse_setting_pre <- c('_', args$type, '.')
collapse_setting <- paste(collapse_setting_pre, collapse="")
plot_name <- paste(name_setting, collapse= collapse_setting)



#dataframe setting
#dataframe selection
data_raw <- read.table(input_file_name, header = T)

if (grepl("\\s", list_file_name)) {data_original_raw <-data_raw
}else{data_list <- read.table(list_file_name, header = F)
    data_raw$count = match(row.names(data_raw), data_list$V1)
    data_original_raw <- subset(na.omit(data_raw), select = -c(count))
    collapse_setting_list_pre <- c(output_file_name, '_gene_selection.txt')
    plot_name_list <- paste(collapse_setting_list_pre, collapse="")
    write.table(data_original_raw, plot_name_list, sep = "\t")
}



#dataframe selection(log and zscore)
data_original_log2 <- log2(data_original_raw+ args$log_number)
data_original_log10 <- log10(data_original_raw+ args$log_number)
data_original_loge <- log(data_original_raw+ args$log_number)
zscore_transformation <- na.omit(t(scale(t(data_original_log2), scale=TRUE, center=TRUE)))
write.table(zscore_transformation, file = "z_score.txt", sep = "\t")
data_original_zscore <- read.table("z_score.txt", header = T)

if (grepl("\\w", args$remove_col[1])){remove <- as.factor(c(args$remove))
    data_raw <-subset(data_raw, select = -c(remove))}

if (args$log == "0"){if (args$zscore == "off"){
			data_original <- data_original_raw
		} else if (args$zscore == "on"){
			data_original <- data_original_zscore
		}
}else if (args$log == "2"){
	data_original <- data_original_log2
}else if (args$log == "10"){
	data_original <- data_original_log10
}else if (args$log == "e"){
	data_original <- data_original_log2
}

data_original2 <-stack(data_original[1:length(colnames(data_original))])
gene <- c(rep(rownames(data_original),length(colnames(data_original))))
data<- data_original2$values
sample<-data_original2$ind
data_set <-data.frame(sample, data, gene)
file.remove("z_score.txt")



#measure time
mid.time <- Sys.time()
mid.time.taken <- mid.time - start.time
paste("Complete dataframe generation! Dataframe generation time (seconds) ; ",round(mid.time.taken, digit = 2))



#axis setting
if (xaxis == "sample"){
    xd<-sample
}else if (xaxis == "gene") {
    xd <- gene
}

if (zaxis == "sample"){
    zd<-sample
}else if (zaxis == "gene") {
    zd <- gene

}

yd <-data



# geom position selection
geom_position <- function() {if (args$geom_position == 1){
                                return ("stack")
                            } else if (args$geom_position == 2){
                                return ("dodge")
                            } else if (args$geom_position == 3){
                                return ("dodge2")
                            } else if (args$geom_position == 4){
                                return ("fill")
                            }
                        }



#color palette selection
#geom color
palette_select_color <- function (){if (args$custom_color[1] != "0"){
                                        return (scale_colour_manual(values =args$custom_color ))
                            }else if (color == 1){
                                        return (scale_colour_hue())
                            }else if (color == 2){
                                return (scale_colour_viridis_d(option = "C"))
                            }else if (color == 3){
                                return (scale_colour_viridis_d(option = "D"))
                            }else if (color == 4){
                                return (scale_colour_grey())
                            }else if (color == 5){
                                return (scale_colour_brewer(palette ="RdBu"))
                            }else if (color == 6){
                                return (scale_colour_brewer(palette ="RdYlBu"))
                            }else if (color == 7){
                                return (scale_colour_brewer(palette ="Reds"))
                            }else if (color == 8){
                                return (scale_colour_brewer(palette ="Blues"))
                            }else if (color == 9){
                                return (scale_colour_brewer(palette ="Paired"))
                            }else if (color == 10){
                                return (scale_colour_brewer(palette = "Set1"))
                            }
                        }

#geom fill
palette_select_fill <- function (){if (args$custom_color[1] != "0"){
                                    return (scale_fill_manual(values =args$custom_color ))
                            }else if (color == 1){
                                return (scale_fill_hue())
                            }else if (color == 2){
                                return (scale_fill_viridis_d(option = "C"))
                            }else if (color == 3){
                                return (scale_fill_viridis_d(option = "D"))
                            }else if (color == 4){
                                return (scale_fill_grey())
                            }else if (color == 5){
                                return (scale_fill_brewer(palette ="RdBu"))
                            }else if (color == 6){
                                return (scale_fill_brewer(palette ="RdYlBu"))
                            }else if (color == 7){
                                return (scale_fill_brewer(palette ="Reds"))
                            }else if (color == 8){
                                return (scale_fill_brewer(palette ="Blues"))
                            }else if (color == 9){
                                return (scale_fill_brewer(palette ="Paired"))
                            }else if (color == 10){
                                return (scale_fill_brewer(palette = "Set1"))
                            }
                        }
#pheatmap color
col_heatmap <- function(){if (args$custom_color[1] != "0"){
                                return (args$custom_color )
                            }else if (color == 1){
                                return (c("blue", "white", "red"))
                            }else if (color == 2){
                                return (c("blue", "light yellow", "red"))
                            }else if (color == 3){
                                return (c("dark green", "white", "red"))
                            }else if (color == 4){
                                return (c("dark green","white" ,"purple"))
                            }else if (color == 5){
                                return (c("black", "white", "red"))
                            }else if (color == 6){
                                return (c("blue", "yellow"))
                            }else if (color == 7){
                                return (c("white", "black"))
                            }else if (color == 8){
                                return (c("white", "red"))
                            }else if (color == 9){
                                return (c("white", "blue"))
                            }else if (color == 10){
                                return (c("white", "green"))
                            }
                }



#ggplot and pheatmap type
#bar_plot
ggplot_bar <- function() {ggplot(data_set, aes(x=xd, y=yd, fill = zd)) +
    geom_col(position = geom_position()) +
    palette_select_fill()+
    labs(x = xaxis, y = "data", fill = zaxis)}

#box_plot
ggplot_box <- function() {ggplot(data_set, aes(x=xd, y=yd, fill = xd)) +
    geom_boxplot(outlier.alpha = 0.3, coef = 1.5) +
    palette_select_fill()+
    labs(x = xaxis, y = "data", fill = xaxis)}

#density_plot
ggplot_density <- function() {ggplot(data_set, aes(x=yd, color = xd)) +
    geom_density(alpha = 1, position = geom_position()) +
    palette_select_color() +
    labs(y = "count", x = "data", color = xaxis)}

#density_fill_plot
ggplot_density_fill <- function() {ggplot(data_set, aes(x=yd, fill = xd)) +
    geom_density(alpha = 0.7, position = geom_position()) +
    palette_select_fill() +
    labs(y = "count", x = "data", fill = xaxis)}

#histogram_plot
ggplot_histogram <- function() {ggplot(data_set, aes(x=yd, fill = xd)) +
    geom_histogram(bin = 30, position = geom_position())+
    palette_select_fill()+
    labs(y = "count", x = "data", fill = xaxis)}

#line_plot
ggplot_line <- function() {ggplot(data_set, aes(x=xd, y = yd, group = zd, color = zd)) +
    geom_line() +
    palette_select_color() +
    labs(x = xaxis, y = "data", color = zaxis)}

#dot_color_plot
ggplot_dot_color <- function() {ggplot(data_set, aes(x=xd, y = yd, color = zd)) +
    geom_point() +
    palette_select_color() +
    labs(x = xaxis, y = "data", color = zaxis)}

#dot_shape_plot
ggplot_dot_shape <- function() {ggplot(data_set, aes(x=xd, y = yd, shape = zd)) +
    geom_point() +
    labs(x = xaxis, y = "data", shape = zaxis)}

#scatter_plot
ggplot_scatter <- function() {ggplot(data_original, aes_string(x= args$scatter_select[1], y = args$scatter_select[2])) +
    geom_point() +
    labs(x = args$scatter_select[1], y = args$scatter_select[2] )}

#violin_plot
ggplot_violin <- function() {ggplot(data_set, aes(x=xd, y=yd, fill = xd)) +
    geom_violin() +
    palette_select_fill()+
    labs(x = xaxis, y = "data", fill = xaxis)}

#heatmap_plot with cluseter on in column and row
pheatmap_heatmap_cc_cr <- function() {pheatmap(data_original, color = colorRampPalette(col_heatmap())(100), fontsize_col= args$letter_size[2], fontsize_row= args$letter_size[2], cluster_rows=T, cluster_cols=T, filename= plot_name)}

#heatmap_plot with cluseter on in column
pheatmap_heatmap_cc <- function() {pheatmap(data_original, color = colorRampPalette(col_heatmap())(100), fontsize_col= args$letter_size[2], fontsize_row= args$letter_size[2], cluster_rows=F, cluster_cols=T, filename = plot_name)}

#heatmap_plot with cluseter on in row
pheatmap_heatmap_cr <- function() {pheatmap(data_original, color = colorRampPalette(col_heatmap())(100), fontsize_col= args$letter_size[2], fontsize_row= args$letter_size[2], cluster_rows=T, cluster_cols=F, filename = plot_name)}

#heatmap_plot with cluseter off
pheatmap_heatmap <- function() {pheatmap(data_original, color = colorRampPalette(col_heatmap())(100), fontsize_col= args$letter_size[2], fontsize_row= args$letter_size[2], cluster_rows=F, cluster_cols=F, filename = plot_name)}



#plot type selection
ggplot_types<- function(){if (args$type == "line"){
    				return (ggplot_line())
			    } else if (args$type == "box"){
                    return (ggplot_box())
                } else if (args$type == "bar"){
                    return (ggplot_bar())#box
                } else if (args$type == "dot_color"){
                    return (ggplot_dot_color())
                } else if (args$type == "dot_shape"){
                    return (ggplot_dot_shape())
                } else if (args$type == "violin"){
                    return (ggplot_violin())
                } else if (args$type == "density"){
                    return (ggplot_density())
                } else if (args$type == "density_fill"){
                    return (ggplot_density_fill())
                } else if (args$type == "histogram"){
                    return (ggplot_histogram())
                } else if (args$type == "scatter"){
                    return (ggplot_scatter())
                } else {
                    return (print("error; graph type(-t, --type) are not provided"))
                }
            }



#limit setting 
limit <- function(){if (args$limit[1] != 0 || args$limit[2] != 0){
                        if (args$type == "scatter" ){lims(x = args$limit, y = args$limit)
                        }else if (args$type == "histogram" || args$type == "density"|| args$type == "density_fill"){lims(x = args$limit)
                        }else {lims(y = args$limit)
                        }
                    }
                }



#axis setting
axis <- function() {if (args$axis_change == "on") {coord_flip()}}



#backgroud selection
back_type <- function(){if (args$style == 1){
                            return (theme_void())#no line with white
                        }else if (args$style == 2){
                            return (theme_classic()) #axis line with white
                        }else if (args$style == 3){
                            return (theme_minimal()) #subline with white
                        }else if (args$style == 4){
                            return (theme_bw()) #grey axis and subline with white
                        }else if (args$style == 5){
                            return (theme_linedraw()) #black axis and subline with white
                        }else if (args$style == 6){
                            return (theme_grey()) #subline with grey
                        }else if (args$style == 7){
                            return (theme_dark()) #subline with black
			}
		}

                            
                            
#legend position and letter size
letter_setting <- theme(legend.position = args$legend_position, axis.text=element_text(size=args$letter_size[1]), axis.title=element_text(size=args$letter_size[2],face="bold"), legend.text=element_text(size=args$letter_size[1]), legend.title = element_text(size=args$letter_size[2],face="bold"))



#figure selection
figure <- function(){if (args$type == "heatmap"){
                            if (args$cluster_select[1] == "on"){
                                if (args$cluster_select[2] == "on"){
                                    pheatmap_heatmap_cc_cr()
                                } else {
                                    pheatmap_heatmap_cc()}
                            }else{
                                if (args$cluster_select[2] == "on"){
                                    pheatmap_heatmap_cr()
                                }else{
                                    pheatmap_heatmap()
                                }
                            }
                    }else {
                        ggplot_types() + limit() + axis() + back_type() + letter_setting + ggsave(plot_name, width =args$plot_size[1], height = args$plot_size[2])
                    }
                }

#start operation
figure()
end.time <- Sys.time()
time.taken <- end.time - start.time
paste("Plotting now! Total operation time (second); ",round(time.taken, digit = 2))
