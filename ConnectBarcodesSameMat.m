function [conn,howManyConnect] = ConnectBarcodesSameMat(out,in,conn,howManyConnect)
whichEdges = [bin2dec(out)+1,bin2dec(in)+1];
conn(sub2ind(size(conn),whichEdges(:,1),whichEdges(:,2))) = 1;
howManyConnect(unique(whichEdges)) = howManyConnect(unique(whichEdges))+1;