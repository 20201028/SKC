Here is the source code of paper *Keyword-Aware Skyline Community Search on Semantics and Structure*

Files extracG.py and extractAtt.py are the code for extracting X%*|V| subgraph
Files large_graph_gen.py and large_attribute_gen.py are the code for generating X * |V| large-scale graphs

### Example of compiling and running (Linux)

Before running, you should:
+ Prepare the dataset and query files according to the templates
+ Check file names and file paths

Then run the following commands,

```shell
$ make
$ main "bk" "0" "3" "l" "3"
# For running queries, five parameters should be specified:
# "bk" is the name of dataset
# "0" means that the similarity threshold
# "3" means that we are running queries on varying experiments
# "l" is the type of attribute, where "i" is the textual attribute, and "l" is the locational attribute 
# "3" is the k value
```






