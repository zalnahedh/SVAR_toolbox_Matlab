function [columnPosition,rowPosition] = excelPositionIrf(info,variablePosition,shockPosition)

startColumn=2;
startRow=2;

columnPosition=startColumn+5*(shockPosition-1);
rowPosition=startRow+(info.T+4)*(variablePosition-1);
end