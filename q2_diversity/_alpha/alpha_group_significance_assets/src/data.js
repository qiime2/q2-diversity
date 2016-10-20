import {
  quantile,
  ascending,
  range,
  min,
  max,
} from 'd3';


function _quartiles(group) {
  return [
    quantile(group, 0.25),
    quantile(group, 0.50),
    quantile(group, 0.75),
  ];
}

function _iqrIndices(group, quartiles) {
  const q1 = quartiles[0];
  const q3 = quartiles[2];
  const k = 1.5;
  const iqr = (q3 - q1) * k;
  let i = 0;
  let j = group.length - 1;
  while (group[i] < (q1 - iqr)) { i += 1; }
  while (group[j] > (q3 + iqr)) { j -= 1; }
  return [i, j];
}

function _toData(group, indices) {
  return indices.map(i => group[i]);
}

export default function setupData(data) {
  const quartiles = [];
  const whiskersIndices = [];
  const whiskersData = [];
  const outliersIndices = [];
  const outliersData = [];
  const minVals = [];
  const maxVals = [];

  // cool naming convention (did i stutter?)
  data.data.data.forEach((vals) => {
    const sortedVals = vals.sort(ascending);

    const quartile = _quartiles(sortedVals);
    quartiles.push(quartile);

    const whiskerIndices = _iqrIndices(sortedVals, quartile);
    whiskersIndices.push(whiskerIndices);

    whiskersData.push(_toData(sortedVals, whiskerIndices));

    const outlierIndices = range(0, whiskerIndices[0])
      .concat(range(whiskerIndices[1] + 1, sortedVals.length));
    outliersIndices.push(outlierIndices);
    outliersData.push(_toData(sortedVals, outlierIndices));

    minVals.push(min(sortedVals));
    maxVals.push(max(sortedVals));
  });
  return {
    quartiles,
    whiskersData,
    outliersData,
    xLabels: data.data.index,
    min: min(minVals),
    max: max(maxVals),
    category: data.category,
    metric: data.metricName,
  };
}
