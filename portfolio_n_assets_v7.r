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

sample_number <- 100

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
rportfolios_tmp <- random_portfolios(port, permutations = 2000, rp_method = "sample")
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

returns.data <- read.table("data_7_assets.txt")
returns.data$VNQI <- NULL
returns.data$VNQ <- NULL
eff.frontier.A <- calculate_frontier(returns.data)
eff.frontier.A.hover_text <- create_hover_text(eff.frontier.A)

returns.data <- read.table("data_7_assets.txt")
eff.frontier.B <- calculate_frontier(returns.data)
eff.frontier.B.hover_text <- create_hover_text(eff.frontier.B)

p <- plot_ly() %>%
  add_trace(x = eff.frontier.A$Risk, y = eff.frontier.A$Return, 
        mode = "markers", type = "scattergl", name = "Portfolio A",
		text = eff.frontier.A.hover_text, hoverinfo = 'text',
        marker = list(size = 5, opacity = 1.0, color = "#000099" )) %>% 
  
  add_trace(x = eff.frontier.B$Risk, y = eff.frontier.B$Return, mode = "markers", 
            type = "scattergl", name = "Portfolio B",
			text = eff.frontier.B.hover_text, hoverinfo = 'text',
            marker = list(color = "#990000", size = 5, opacity = 1.0)) %>% 
			
  layout(title = "Random Portfolios with Plotly",
         yaxis = list(title = "Mean Returns", tickformat = ".2%"),
         xaxis = list(title = "Standard Deviation", tickformat = ".2%"),
         plot_bgcolor = "#FFFFFF",
         paper_bgcolor = "#FFFFFF"
         )
p
