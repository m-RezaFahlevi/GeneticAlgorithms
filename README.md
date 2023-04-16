# Genetic Algorithm
___
[![Licence](https://img.shields.io/github/license/Ileriayo/markdown-badges?style=for-the-badge)](./LICENSE)

The main parts of this repository is in folder GeneticAlgorithm. SteadyStateGeneticAlgorithm.R and steady_state_genetic_algorithm.py are the main source code. **Real-Parameter Optimazation Problem** are supposed to be solved in this repository, at least for now.

<p align="center">
<img alt="arg min spherefun"
     src="https://github.com/m-RezaFahlevi/GeneticAlgorithms/blob/main/GeneticAlgorithm/visualization/obt_data1681577275.gif"/>
     <caption><br/>Solving $\arg \min_{x_1,x_2\in\mathbb{R}^2}\phi(x_1, x_2) = x_1^2 + x_2^2$</caption>
</p>

## SteadyStateGeneticAlgorithm.R

For file SteadyStateGeneticAlgorithm.R, there are 2 library must be installed first, **ggplot2 and dplyr**. These library can be installed via RStudio or execute following code in your R's terminal.

```{r}
install.packages(c("ggplot2", "dplyr"))
```

When the execution of file SteadyStateGeneticAlgorithm.R is finish, it's show that max[f(x,y)] if and only if x,y = 0.

## steady_state_genetic_algorithm.py

![](https://github.com/m-RezaFahlevi/GeneticAlgorithms/blob/main/GeneticAlgorithm/www/Screenshot-20210623154036-1365x735.png)

### Prelude
Install these 2 python module, **numpy** and **pandas**
```
pip install numpy
```
```
pip install pandas
```
### How to Execute It
Run the file in your python IDE, or, if python had been installed in your system, open your terminal, change the current working directory to the directory where file steady_state_genetic_algorithm.py belongs to, then execute the following code

```
python steady_state_genetic_algorithm.py
```
