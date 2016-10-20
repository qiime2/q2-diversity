export function setupXLabel(svg, width, height, label, xAxis) {
  let maxHeight = 0;

  svg.select('.x.axis')
    .attr('transform', `translate(0,${height})`)
    .call(xAxis)
    .selectAll('text')
    .style('text-anchor', 'end')
    .attr('dx', '-.8em')
    .attr('dy', '-0.5em')
    .attr('transform', function setHeight() {
      const textHeight = this.getComputedTextLength();
      if (textHeight > maxHeight) maxHeight = textHeight;
      return 'rotate(-90)';
    });

  svg.select('.x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(${(width / 2)},${(height + maxHeight + 40)})`)
    .text(label);

  return maxHeight;
}

export function setupYLabel(svg, height, label, yAxis) {
  svg.select('.y.axis')
    .call(yAxis);

  svg.select('.y.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(-50,${(height / 2)})rotate(-90)`)
    .text(label);
}
