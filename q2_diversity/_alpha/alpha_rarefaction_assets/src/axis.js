import { max } from 'd3';

export function setupXLabel(svg, width, height, label, xAxis) {
  svg.select('#chart .x.axis')
    .attr('transform', `translate(0,${height})`)
    .transition()
    .call(xAxis);

  svg.select('#chart .x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .text(label)
    .attr('transform', `translate(${(width / 2)},${(height + 30)})`);

  svg.select('#subChart .x.axis')
    .attr('transform', `translate(0,${height})`)
    .transition()
    .call(xAxis);

  svg.select('#subChart .x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .text(label)
    .attr('transform', `translate(${(width / 2)},${(height + 30)})`);
}

export function setupYLabels(svg, height, label, yAxisChart, yAxisSubChart) {
  // for some reason using transition here breaks the text calculation below
  function setupLabel(selection, l) {
    const setup = selection.attr('text-anchor', 'middle')
      .style('font', '12px sans-serif')
      .text(l);
    return setup;
  }
  const chartLine = svg.select('#chart .y.axis')
    .call(yAxisChart);

  const chartLabel = svg.select('#chart .y.label').call(setupLabel, label);
  const all = Array.from(chartLine.selectAll('text')._groups[0]).map(d => d.getComputedTextLength());
  const textHeight = max(all) + 20;
  chartLabel.attr('transform', `translate(-${textHeight},${(height / 2)})rotate(-90)`);

  svg.select('#subChart .y.axis')
    .call(yAxisSubChart);
  const subChartLabel = svg.select('#subChart .y.label').call(setupLabel, 'samples');
  subChartLabel.attr('transform', `translate(-${textHeight},${(height / 2)})rotate(-90)`);

  return textHeight;
}
