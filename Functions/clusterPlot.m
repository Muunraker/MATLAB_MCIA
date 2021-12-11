function x = clusterPlot(global_xcoords,global_ycoords,clusList)
% Plots a simple scatterplot given a list of coordinates and clusters with unique symbols
% Clusters must be given in form of a numeric array, e.g. [1, 3, 1, 1, 3, 3, 4,...]
% Inputs: 
%         * global_xcoords    list of x-coordinates
%         * global_ycoords    list of y-coordinates
%         * clusList          numeric array of corresponding clusters for each coordinate

    x = [];

    markerList = {'*','.','o','x','+'};
    colorList = {'b','r','g','m','y','c'};
    clusNames = unique(clusList);
    numClusters = length(clusNames);
    
    if numClusters > length(markerList)*length(colorList)
       error('Number of clusters exceeds number of plotting markers') 
    end

    for i=0:numClusters-1
        clusXcoords = global_xcoords(clusList==clusNames(i+1));
        clusYcoords = global_ycoords(clusList==clusNames(i+1));

        scatter(clusXcoords,clusYcoords,[colorList{1+mod(i,length(colorList))},...
        markerList{1+floor(i/length(colorList))}]); hold on
    end
end

