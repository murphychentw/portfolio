# read assets from file
# read portfolio from file
# plot return vs risk chart
# plot yearly return chart
# plot cumulative return chart

library(PortfolioAnalytics)
library(quantmod)
library(PerformanceAnalytics)
library(zoo)
library(plotly)
library(ROI)
library(DEoptim)
library(doParallel)
require(ROI.plugin.glpk)
require(ROI.plugin.quadprog)

registerDoParallel()

sample_number <- 200

calculate_frontier <- function(returns.data) {

  assets_number <- ncol(returns.data)
  assets_name = colnames(returns.data)

  returns.data <- na.omit(returns.data)

  # Save mean return vector and sample covariance matrix
  meanReturns <- colMeans(returns.data)
  covMat <- cov(returns.data)

  # Start with the names of the assets
  port <- portfolio.spec(assets = assets_number, category_labels = assets_name)

  # Box
  port <- add.constraint(port, type = "box", min = 0.01, max = 1.0)
 
  # Leverage
  port <- add.constraint(portfolio = port, type = "weight_sum", min_sum=0.99, max_sum=1.01)

  # Generate random portfolios
  rportfolios_tmp <- random_portfolios(port, permutations = 3000, rp_method = "sample")
  rportfolios = rportfolios_tmp[1:sample_number,]

  # Get minimum variance portfolio
  minvar.port <- add.objective(port, type = "risk", name = "var")

  # Optimize
  minvar.opt <- optimize.portfolio(returns.data, minvar.port, optimize_method = "random", rp = rportfolios)
								 
  # Generate maximum return portfolio
  maxret.port <- add.objective(port, type = "return", name = "mean")

  # Optimize
  maxret.opt <- optimize.portfolio(returns.data, maxret.port, optimize_method = "random", rp = rportfolios)

  # Generate vector of returns
  minret <- 0.06/100
  maxret <- maxret.opt$weights %*% meanReturns
 
  vec <- seq(minret, maxret[1], length = sample_number)

  eff.frontier <- data.frame(Risk = rep(NA, sample_number),
    Return = rep(NA, sample_number))
						   
  frontier.weights <- mat.or.vec(nr = sample_number, nc = ncol(returns.data))
  colnames(frontier.weights) <- colnames(returns.data)

  for(i in 1:sample_number){
    eff.port <- port
    eff.port <- add.constraint(eff.port, type = "return", name = "mean", return_target = vec[i])
    eff.port <- add.objective(eff.port, type = "risk", name = "var")
    eff.opt <- optimize.portfolio(returns.data, eff.port, optimize_method = "random", rp = rportfolios)
  
    eff.frontier$Risk[i] <- sqrt(t(eff.opt$weights) %*% covMat %*% eff.opt$weights)
    eff.frontier$Return[i] <- eff.opt$weights %*% meanReturns
    frontier.weights[i,] = eff.opt$weights
  
    print(paste(round(i/sample_number * 100, 0), "% done..."))
  }

  cbind(eff.frontier, frontier.weights)
}

create_hover_text <- function(eff.frontier) {
  hover_text <- array(NA, c(sample_number))
  for(i in 1:sample_number){
    hover_text[i] <- ""
    for (j in 1:ncol(eff.frontier)) {
      hover_text[i] <- paste(hover_text[i], colnames(eff.frontier[j]), ": ", eff.frontier[i,j])
	  hover_text[i] <- paste(hover_text[i], "", sep="\n")
    }
  }
  hover_text
}

######## plot portfolio return vs risk chart ########

returns.data <- read.table("data_8_assets.txt")
returns.data$VNQI <- NULL
returns.data$VNQ <- NULL
eff.frontier.A <- calculate_frontier(returns.data)
eff.frontier.A.hover_text <- create_hover_text(eff.frontier.A)

returns.data <- read.table("data_8_assets.txt")
eff.frontier.B <- calculate_frontier(returns.data)
eff.frontier.B.hover_text <- create_hover_text(eff.frontier.B)

meanReturns <- colMeans(returns.data)
covMat <- cov(returns.data)

portfolio <- read.table("portfolio_8_assets.txt")
portfolio_number <- nrow(portfolio)

port <- data.frame(Risk = rep(NA, portfolio_number), Return = rep(NA, portfolio_number))
for(i in 1:portfolio_number){
  weights <- t(portfolio[i,])
  port$Return[i] <- meanReturns %*% weights
  port$Risk[i] <-  sqrt(t(weights) %*% covMat %*% weights)
  #port$Risk[i] <- sd(t(t(weights) %*% t(returns.data)))
}

p <- plot_ly() %>%
  layout(title = "Portfolio Optimization with R and Plotly",
    yaxis = list(title = "Mean Returns", tickformat = ".2%"),
    xaxis = list(title = "Standard Deviation", tickformat = ".2%"),
    plot_bgcolor = "#FFFFFF",
    paper_bgcolor = "#FFFFFF")

p <- add_trace(p, x = eff.frontier.A$Risk, y = eff.frontier.A$Return, 
    mode = "markers", type = "scattergl", name = "Annuity+BNDX+BND+\nVXUS+VTI",
	text = eff.frontier.A.hover_text, hoverinfo = 'text', marker = list(size = 5))
  
p <- add_trace(p, x = eff.frontier.B$Risk, y = eff.frontier.B$Return, mode = "markers", 
    type = "scattergl", name = "Annuity+BNDX+BND+\nVXUS+VTI+VNQI+VNQ",
	text = eff.frontier.B.hover_text, hoverinfo = 'text', marker = list(size = 5))

p <- add_trace(p, x = port$Risk, y = port$Return,
    mode = "markers", type = "scattergl", name = rownames(portfolio), marker = list(size = 8))

#p.ay <- floor(runif(portfolio_number, +2, +6))*10
p.ay <- array(1:5, dim=c(1, portfolio_number))*15

p <- add_annotations(p, x = port$Risk, y = port$Return,
    text = rownames(portfolio), xanchor = "left", font=list(size = 10), ax = p.ay,
	ay = p.ay, arrowhead = 2, arrowsize = 1, arrowwidth = 1)
	
port.risk.return <- p

######## plot portfolio yearly return chart (line chart) ########

returns.data <- as.matrix(returns.data)
x.axis <- factor(rownames(returns.data), levels = rownames(returns.data))
data <- data.frame(x.axis)
data$x.axis <- factor(data$x.axis, levels = data[["x.axis"]])
pure_asset_number <- 8

p <- plot_ly(data) %>%
  layout(title = "Portfolio Yearly Return",
    yaxis = list(title = "Return", tickformat = ".2%"), xaxis = list(title = "Year"),
    plot_bgcolor = "#FFFFFF", paper_bgcolor = "#FFFFFF")

# pure asset
for(i in 1:pure_asset_number){
  p <- add_trace(p, x = ~x.axis, y = returns.data %*% t(portfolio[i,]), 
    type = "scatter", mode = "lines", name = rownames(portfolio)[i],
	line = list(width = 0.5, dash = "dot"))
}

# portfolios
for(i in (pure_asset_number+1):portfolio_number){
  p <- add_trace(p, x = ~x.axis, y = returns.data %*% t(portfolio[i,]), 
    type = "scatter", mode = "lines", name = rownames(portfolio)[i],
    line = list(width = 1))
}

port.yearly.return.line.chart <- p

# text version:
write.table(returns.data %*% t(portfolio[]), "port.yearly.return.txt", sep="\t")

######## plot portfolio cumulative return chart (line chart) ########

p <- plot_ly(data) %>%
  layout(title = "Portfolio Cumulative Return",
    yaxis = list(title = "Return", tickformat = ".2%"), xaxis = list(title = "Year"),
    plot_bgcolor = "#FFFFFF", paper_bgcolor = "#FFFFFF")

cumulative.return <- returns.data %*% t(portfolio[])
	
# pure asset
number_of_years <- nrow(returns.data)
for(i in 1:pure_asset_number){
  weights <- t(portfolio[i,])
  cumulative.return[1,i] <- returns.data[1,] %*% weights
  for (j in 2:number_of_years) {
    cumulative.return[j,i] <- (cumulative.return[j-1,i]+1)*((returns.data[j,] %*% weights)+1)-1
  }
  p <- add_trace(p, x = ~x.axis, y = cumulative.return[,i], 
    type = "scatter", mode = "lines", name = rownames(portfolio)[i],
	line = list(width = 0.5, dash = "dot"))
}

# portfolios
for(i in (pure_asset_number+1):portfolio_number){
  weights <- t(portfolio[i,])
  cumulative.return[1,i] <- returns.data[1,] %*% weights
  for (j in 2:number_of_years) {
    cumulative.return[j,i] <- (cumulative.return[j-1,i]+1)*((returns.data[j,] %*% weights)+1)-1
  }
  p <- add_trace(p, x = ~x.axis, y = cumulative.return[,i], 
    type = "scatter", mode = "lines", name = rownames(portfolio)[i],
	line = list(width = 1))
}

port.cumulative.return.line.chart <- p

# text version:
write.table(cumulative.return, "port.cumulative.return.txt", sep="\t")

######## covariance matrix ########

write.table(cov(returns.data), "port.covariance.txt", sep="\t")

######## corelation matrix ########

write.table(cor(returns.data), "port.corelation.txt", sep="\t")

######## histogram of return ########

port.returns <- returns.data %*% t(portfolio[])

p <- plot_ly(data, alpha = 0.6) %>%
  layout(title = "Histogram of Return",
    yaxis = list(title = "Times"), xaxis = list(title = "Return"),
    plot_bgcolor = "#FFFFFF", paper_bgcolor = "#FFFFFF", barmode = "stack")

for(i in 1:pure_asset_number){
  p <- add_trace(p, x = port.returns[,i], type = "histogram", name = rownames(portfolio)[i])
}
	
for(i in (pure_asset_number+1):portfolio_number){
  p <- add_trace(p, x = port.returns[,i], type = "histogram", name = rownames(portfolio)[i])
}
