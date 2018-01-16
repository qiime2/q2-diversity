import {
  select,
  scaleBand,
  scaleLinear,
  axisBottom,
  axisLeft,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';
import plotBoxes from './bar';


export function render(svg, data) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 70, right: 50, bottom: 50 };
  const chart = svg.select('g');

  const { xLabels, min, max, column, metric } = data;

  const x = scaleBand().padding(0.1).domain(xLabels).range([0, width]);
  const y = scaleLinear().domain([min, max]).range([height, 0]).nice();

  const xAxis = axisBottom();
  const yAxis = axisLeft();

  xAxis.scale(x);
  yAxis.scale(y);

  chart.attr('transform', `translate(${margin.left},${margin.top})`);

  const maxLabelX = setupXLabel(svg, width, height, column, xAxis);
  setupYLabel(svg, height, metric, yAxis);

  plotBoxes(chart, data, x, y);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top + maxLabelX);
}

export function kwStats(body, data) {
  select('#kw-all-h').text(data.kwAll.H);
  select('#kw-all-p').text(data.kwAll.p);
  select('#kw-pairwise').html(data.kwPairwise);
  select('#kw-csv').attr('href', data.kwCSVPath);
}

export function warnings(body, data) {
  const { filtered: { initial, filtered } } = data;
  if (initial !== filtered) {
    select('#filtered-samples')
      .style('display', null)
      .html(`Some samples were filtered from the input alpha diversity data
        because they were missing metadata values.<strong>The input alpha
        diversity data contained ${initial} samples, but the analysis was run
        on only ${filtered} samples.</strong>`);
  } else {
    select('#filtered-samples').style('display', 'none').html(null);
  }
}
