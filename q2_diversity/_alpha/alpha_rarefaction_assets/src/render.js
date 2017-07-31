import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory10,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';

function renderPlot(svg, data, x, y, category) {
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

  chart.selectAll('.legend').remove();
  const legend = chart.append('g')
    .attr('class', 'legend')
    .attr('x', 0)
    .attr('y', 0)
    .attr('overflow', 'scroll')
    .style('border', '3px solid black')
    .attr('height', 100)
    .attr('width', 110);
  legend.selectAll('g')
        .data(Array.from(setGroups))
      .enter()
        .append('text')
          .attr('x', 25)
          .attr('y', (d, i) => (i * 25) + 8)
          .attr('height', 15)
          .attr('width', 100)
          .style('fill', d => color(d))
          .text(d => d);
}

export default function render(svg, data, category) {
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

  renderPlot(svg, data, x, y, category);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}

// export function stats(body, data) {
//   const { stats: { method, testStat, pVal, sampleSize } } = data;
//   select('#method').text(method);
//   select('#test-stat').text(testStat);
//   select('#p-val').html(pVal);
//   select('#sample-size').html(sampleSize);
// }
