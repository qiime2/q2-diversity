export function setupXLabel(svg, width, height, label, xAxis) {
  svg.select('.x.axis')
    .attr('transform', `translate(0,${height})`)
    .call(xAxis);

  const l = svg.select('.x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(${(width / 2)},${(height + 40)})`)
    .text(label);
  const textHeight = l.node().getComputedTextLength();
  return Math.max(textHeight, 30);
}

export function setupYLabel(svg, height, label, yAxis) {
  svg.select('.y.axis')
    .call(yAxis);

  const l = svg.select('.y.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(-50,${(height / 2)})rotate(-90)`)
    .text(label);
  const textHeight = l.node().getComputedTextLength();
  return Math.max(textHeight, 30);
}
