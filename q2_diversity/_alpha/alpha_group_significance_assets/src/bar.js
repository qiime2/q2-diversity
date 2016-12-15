export default function plotBoxes(chart, data, x, y) {
  const { quartiles, outliersData, whiskersData, xLabels } = data;
  const halfWidth = x.bandwidth() / 2;
  const quarterWidth = halfWidth / 2;

  const centerUpdate = chart.selectAll('line.center').data(whiskersData, d => d);
  centerUpdate.exit().remove();
  const centerEnter = centerUpdate.enter().append('line');
  centerUpdate.merge(centerEnter)
    .attr('class', 'center')
    .attr('x1', (_, i) => (x(xLabels[i]) + halfWidth))
    .attr('y1', d => y(d[0]))
    .attr('x2', (_, i) => (x(xLabels[i]) + halfWidth))
    .attr('y2', d => y(d[1]))
    .attr('stroke-dasharray', '2,2')
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const outlierUpdate = chart.selectAll('.outlier').data(outliersData, d => d);
  outlierUpdate.exit().remove();
  const outlierEnter = outlierUpdate.enter().append('g').attr('class', 'outlier');
  const outliers = outlierUpdate.merge(outlierEnter)
    .attr('transform', (_, i) => (`translate(${x(xLabels[i]) + halfWidth},0)`));

  const circleUpdate = outliers.selectAll('circle').data(d => d);
  circleUpdate.exit().remove();
  const circleEnter = circleUpdate.enter().append('circle');
  circleUpdate.merge(circleEnter)
    .attr('r', 4)
    .attr('stroke-width', 1)
    .attr('stroke', 'black')
    .attr('stroke-opacity', 0.3)
    .attr('fill', 'black')
    .attr('fill-opacity', 0.33)
    .attr('cx', 0)
    .attr('cy', d => y(d));

  const boxUpdate = chart.selectAll('rect.box').data(quartiles, d => d);
  boxUpdate.exit().remove();
  const boxEnter = boxUpdate.enter().append('rect');
  boxUpdate.merge(boxEnter)
    .attr('class', 'box')
    .attr('x', (_, i) => x(xLabels[i]))
    .attr('y', d => y(d[2]))
    .attr('width', x.bandwidth())
    .attr('height', d => (y(d[0]) - y(d[2])))
    .attr('fill', 'white')
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const medianUpdate = chart.selectAll('line.median').data(quartiles, d => d);
  medianUpdate.exit().remove();
  const medianEnter = medianUpdate.enter().append('line');
  medianUpdate.merge(medianEnter)
    .attr('class', 'median')
    .attr('x1', (_, i) => x(xLabels[i]))
    .attr('y1', d => y(d[1]))
    .attr('x2', (_, i) => (x(xLabels[i]) + x.bandwidth()))
    .attr('y2', d => y(d[1]))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const whiskerUpdate = chart.selectAll('.whisker').data(whiskersData, d => d);
  whiskerUpdate.exit().remove();
  const whiskerEnter = whiskerUpdate.enter().append('g').attr('class', 'whisker');
  const whiskers = whiskerUpdate.merge(whiskerEnter)
    .attr('transform', (_, i) => (`translate(${x(xLabels[i]) + quarterWidth},0)`));

  const whiskerLines = whiskers.selectAll('line').data(d => d);
  whiskerLines.enter().append('line')
    .attr('class', 'median')
    .attr('x1', 0)
    .attr('y1', d => y(d))
    .attr('x2', halfWidth)
    .attr('y2', d => y(d))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');
}
