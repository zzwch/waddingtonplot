# Waddington Epigenetic Landscape
> Waddington’s epigenetic landscape is probably the most famous and most powerful metaphor in developmental biology. Cells, represented by balls, roll downhill through a landscape of bifurcating valleys. Each new valley represents a possible cell fate and the ridges between the valleys maintain the cell fate once it has been chosen. (Quoted from [reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3372930/))

# [R package waddingtonplot](https://github.com/zzwch/waddingtonplot) 
Plot Waddington Epigenetic Landscape using ggridges, inspired by Figure 1 in [reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3372930/).

<img src= "https://raw.githubusercontent.com/lizc07/myScripts/master/images/waddington.ref.png" width = 400> vs <img src="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3372930/bin/nihms376417f1.jpg" width = 400>

# Install and try it
```r
devtools::install_github("zzwch/waddingtonplot")
waddingtonplot::waddingtonplot()
```

# custom your graph 
add external images using 
```r
waddingtonplot::waddingtonplot(c(1,2,4,8), do.return = T) +
annotation_custom(grob = grid::rasterGrob(image = EBImage::readImage(files = "img_20180217191752.jpg")),
                    450, 550, 77, 85)
```
<img src= "https://raw.githubusercontent.com/lizc07/myScripts/master/images/waddington.toy.png" width = 600>
