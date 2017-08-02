import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory10,
  select,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';

function renderPlot(svg, data, x, y, category, legend) {
  const chart = svg.select('g');

  const depthIndex = data.data.columns.indexOf('depth');
  const medianIndex = data.data.columns.indexOf('median');
  let groupIndex = data.data.columns.indexOf('sample-id');
  if (groupIndex === -1) {
    groupIndex = data.data.columns.indexOf(category);
  }
  const points = [data.data.data.sort((a, b) => a[depthIndex] - b[depthIndex])][0];
  const setGroups = new Set(Array.from(points, d => d[groupIndex]));
  const color = scaleOrdinal(schemeCategory10)
    .domain(setGroups);

  chart.selectAll('circle').remove();
  chart.selectAll('dot')
      .data(points)
    .enter()
      .append('circle')
        .attr('cx', d => x(d[depthIndex]))
        .attr('cy', d => y(d[medianIndex]))
        .attr('r', 4)
        .style('stroke', d => color(d[groupIndex]))
        .style('fill', d => color(d[groupIndex]));

  legend.selectAll('.legend').remove();
  const arrGroups = Array.from(setGroups);
  legend.attr('height', arrGroups.length * 20);
  let ly = 0;
  for (const [i, entry] of arrGroups.entries()) {
    ly = i * 20;
    const c = color(entry);
    legend.append('rect')
      .attr('class', 'legend')
      .attr('x', 0)
      .attr('y', ly - 2.5)
      .attr('width', 15)
      .attr('height', 5)
      .style('stroke', 'darkGrey')
      .style('fill', 'white');
    legend.append('circle')
      .attr('class', 'legend')
      .attr('cx', 30)
      .attr('cy', ly)
      .attr('r', 5)
      .style('stroke', 'darkGrey')
      .style('fill', c);
    legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .style('font', '10px sans-serif')
      .text(entry);
  }
  select(legend.node().parentNode).attr('viewBox', `0 0 200 ${ly + 10}`);
}

export default function render(svg, data, category, legend) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 70, right: 50, bottom: 50 };
  const chart = svg.select('g');

  const { xAxisLabel, yAxisLabel, minX, maxX, minY, maxY } = data;

  const xAxis = axisBottom();
  const yAxis = axisLeft();

  let pad = (maxX - minX) * 0.03;
  if (Number.isInteger(minX) && Number.isInteger(maxX)) {
    pad = Math.max(Math.round(pad), 1);
    const between = Math.max(3, (maxX - minX) + (2 * pad));
    xAxis.ticks(Math.min(between, 12), 'd');
  }

  const x = scaleLinear().domain([minX - pad, maxX + pad]).range([0, width]).nice();
  const y = scaleLinear().domain([minY, maxY]).range([height, 0]).nice();

  xAxis.scale(x);
  yAxis.scale(y);

  chart.attr('transform', `translate(${margin.left},${margin.top})`);

  setupXLabel(svg, width, height, xAxisLabel, xAxis);
  setupYLabel(svg, height, yAxisLabel, yAxis);

  renderPlot(svg, data, x, y, category, legend);

  svg.attr('width', 1400 + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}
