#!/usr/bin/env Rscript

library(ggplot2)

load("../output/article-data-1.rda")

print(group_resolution_data)

q1 <- function(x) {quantile(x, 1/4)}
q3 <- function(x) {quantile(x, 3/4)}
se <- function(x) {sd(x) / sqrt(length(x))}
ci <- function(x) {se(x) * qt(0.975, length(x) - 1)}

generate_plot <- function(dat, direction, same_group) {
	direction_diff <- paste0(direction, "_diff")
	dat <- dplyr::filter(dat, direction == !!direction, same_group == !!same_group)
	dat_summary <- dat |>
		dplyr::group_by(.data[[direction_diff]]) |>
		dplyr::summarise(mean = mean(value), upper = mean(value) + ci(value), lower = mean(value) - ci(value))
	
	ggplot() +
		geom_ribbon(data = dat_summary, aes(x = .data[[direction_diff]], ymin = upper, ymax = lower), fill = "grey") +
		geom_line(data = dat_summary, aes(x = .data[[direction_diff]], y = mean)) +
		geom_point(data = dat, aes(x = .data[[direction_diff]], y = value))
}

generate_plot(group_resolution_data_long, direction = "row", same_group = FALSE)
generate_plot(group_resolution_data_long, direction = "col", same_group = FALSE)
generate_plot(group_resolution_data_long, direction = "row", same_group = TRUE)
generate_plot(group_resolution_data_long, direction = "col", same_group = TRUE)

row_diff_data <- group_resolution_data |>
	dplyr::group_by(row_diff) |>
	dplyr::summarise_all(list(mean = mean, q1 = q1, q3 = q3))

col_diff_data <- group_resolution_data |>
	dplyr::group_by(col_diff) |>
	dplyr::summarise_all(list(mean = mean, q1 = q1, q3 = q3))

ggplot() +
	geom_ribbon(data = row_diff_data, aes(x = row_diff, ymin = row_1_3_q1, ymax = row_1_3_q3), fill = "grey") +
	geom_line(data = row_diff_data, aes(x = row_diff, y = row_1_3_mean)) +
	geom_point(data = group_resolution_data, aes(x = row_diff, y = row_1_3, col = col_diff))

ggplot() +
	geom_ribbon(data = col_diff_data, aes(x = col_diff, ymin = col_1_3_q1, ymax = col_1_3_q3), fill = "grey") +
	geom_line(data = col_diff_data, aes(x = col_diff, y = col_1_3_mean)) +
	geom_point(data = group_resolution_data, aes(x = col_diff, y = col_1_3, col = row_diff))

ggplot() +
	geom_ribbon(data = row_diff_data, aes(x = row_diff, ymin = row_1_2_q1, ymax = row_1_2_q3), fill = "grey") +
	geom_line(data = row_diff_data, aes(x = row_diff, y = row_1_2_mean)) +
	geom_point(data = group_resolution_data, aes(x = row_diff, y = row_1_2, col = col_diff))

ggplot() +
	geom_ribbon(data = col_diff_data, aes(x = col_diff, ymin = col_1_2_q1, ymax = col_1_2_q3), fill = "grey") +
	geom_line(data = col_diff_data, aes(x = col_diff, y = col_1_2_mean)) +
	geom_point(data = group_resolution_data, aes(x = col_diff, y = col_1_2, col = row_diff))

model <- lm(row_1_3 ~ col_diff * row_diff, data = group_resolution_data)
summary(model)
anova(model)
