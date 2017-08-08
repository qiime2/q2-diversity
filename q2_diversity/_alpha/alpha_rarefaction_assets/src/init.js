/* global categories */
/* global metrics */
import { select } from 'd3';
import state from './state';
import { addMetricPicker, addCategoryPicker } from './toolbar';

export default function init() {
  const metric = metrics[0];
  const category = categories[0];
  // body
  const body = select('#main');
  // controls row, plot row. Order matters here.
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const controlsRow = body.insert('div', ':first-child').attr('class', 'viz row');
  // within controls row we have controls
  const controlsDiv = controlsRow.append('div').attr('class', 'col-lg-12');
  // within plot row we have plot and legend
  const plotCol = plotRow.append('div')
                    .attr('class', 'col-lg-9')
                    .style('width', '1120px');
  const legendCol = plotRow.append('div')
                      .style('height', '470px')
                      .attr('class', 'col-lg-2');
  // within plot col we have plot svg
  const plotSvg = plotCol.append('svg')
                    .attr('viewBox', '0 0 1120 470');
  // within plot svg we have chart g
  const chart = plotSvg.append('g');
  // within chart g we have x axis, y axis, & axis labels
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');
  // within legendCol we have legend title and legend box
  const legendTitle = legendCol.append('div')
                        .style('height', '25px')
                        .style('width', '300px')
                        .style('overflow-y', 'hidden')
                        .style('overflow-x', 'auto')
                        .append('svg')
                          .attr('viewBox', '0 0 200 20')
                          .attr('class', 'legendTitleSvg')
                          .append('g')
                            .attr('class', 'legendTitle');
  const legendBox = legendCol.append('div')
                    .style('height', '445px')
                    .style('width', '300px')
                    .style('overflow-y', 'scroll')
                    .style('overflow-x', 'auto')
                    .attr('class', 'legendBoxDiv')
                    .append('svg')
                      .attr('viewBox', '0 0 300 445')
                      .attr('class', 'legendBoxSvg')
                      .append('g');
  body.insert('h1', ':first-child').text('Alpha Rarefaction');
  // D3
  state.initialize(metric, category, controlsDiv, plotSvg, legendBox, legendTitle);
  addMetricPicker(controlsDiv, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsDiv, categories, category);
  }
}
