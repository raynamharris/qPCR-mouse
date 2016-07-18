#install.packages("gapminder") # you only need to do this once
library(gapminder)
gapminder = gapminder # make an object called "gapminder". you could call this whatever you want
head(gapminder) # look at the first lines
str(gapminder) # look at the structure
names(gapminder) # look at the names of the columns
head(gapminder[1]) # look at the first bit of the first column
summary(gapminder) # summarize the data

# plot the data
plot(lifeExp ~ gdpPercap, data = gapminder, col = continent, main = "My Cool Plot", ylab = "Life Exp", xlab = "GDP")

# another way to plot the data. make sure you put x and y in the right order. This way expects x,y. Using the tilde expects y ~ x.
#plot(gapminder$lifeExp, gapminder$gdpPercap)

# how to get help
?plot

library(ggplot2)

# plot, color by continent
ggplot(data = gapminder, 
       aes(x = year, y = lifeExp, color = continent))+
       geom_point()+
       geom_line()

# plot, color only the lines by continent
ggplot(data = gapminder, 
       aes(x = year, y = lifeExp))+
      geom_point()+
      geom_line(aes(color = continent))+
      theme_bw()

# plot, color only the points by continent
ggplot(data = gapminder, 
       aes(x = year, y = lifeExp))+
  geom_point(aes(color = continent))+
  geom_line()+
  theme_bw()

# plot each country as its own bar graph

# subset the data for only africa
gapminder_africa = gapminder[gapminder$continent=="Africa",] # this says "make an object called gapminder_africa that contains only the rows of gapminder where the continent is africa"
# plot by country within africa
ggplot(data = gapminder_africa, 
    aes(x = year, y = lifeExp))+
    geom_line()+
    theme_bw()+
    facet_wrap(~country)

# using dplyr to filter data frames
library(dplyr)
summary(gapminder)
africa = filter(gapminder,continent=="Africa") # this says "make an object called africa that contains only the rows of gapminder where the continent is africa"
head(africa)

# plot by country within africa
ggplot(data = africa, 
       aes(x = year, y = lifeExp))+
  geom_line()+
  theme_bw()+
  facet_wrap(~country)

europe = filter(gapminder, continent=="Europe")

# plot this to save using Quartz
fig1 = ggplot(data = europe, 
       aes(x = year, y = lifeExp))+
  geom_line()+
  theme_bw()+
  facet_wrap(~country)

quartz() # on macs, you need to say "windows()" on PCs
fig1
