# Multiple-Algorithm Coastal wetland Eco-geomorphology Simulator

**Multiple-Algorithm Coastal wetland Eco-geomorphology Simulator (MACES)** is a model framework to assess the impact of structural uncertainty of eco-geomorphology models on the prediction of coastal wetland evolution under intensified natural and human-induced disturbances. The MACES model consists of two components: a one-dimensional (1-D) transect-based hydrodynamics module (MACES-hydro) and four eco-geomorphology modules with multiple algorithm implementations (MACES-geomor). MACES-hydro simulates water level, tide velocity, significant wave height, bottom shear stress, suspended sediment and other hydrodynamic conditions in a coastal transect along the elevation and land cover gradient with water level and wind speed at the seaward side as inputs. Based on the simulated hydrodynamics, MACES-geomor calculates sediment deposition and OM burial at each grid cell of the coastal transect and lateral erosion at the wetland edge grid cell. At the end of each year, MACES updates the transect elevation profile and land cover. A basic feature of MACES is that different combinations of algorithms within four eco-geomorphology modules can be configured to test different model structures and evaluate their performances.


You can use the [editor on GitHub](https://github.com/tanzeli1982/MACES/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/tanzeli1982/MACES/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
